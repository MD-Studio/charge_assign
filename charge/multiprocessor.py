import multiprocessing as mp
import os
import sys
from typing import Any, Generator, List, Tuple, Type

from charge.util import print_progress

class _Stop:
    pass


class MultiProcessor:
    """
    A helper for processing data using multiple processes.

    This somewhat resembles multiprocessing.Pool, but works \
    differently. To use it, you write a class with a process() method. \
    You pass this class to the MultiProcessor when you construct the \
    MultiProcessor, together with a tuple of arguments to pass to your \
    class's constructor.

    MultiProcessor will then start a number of processes, and create \
    an instance of your class in each process. It will then wait.

    When you then call processed() with a list of tuples of arguments, \
    MultiProcessor will, for each tuple in the list, call process() on \
    one of the instances of your class, with the contents of the tuple \
    as arguments, and yield the value it returns back to you. Note that \
    because the data are processed in parallel, the results may be out \
    of order!

    You can call processed() multiple times to process multiple lists \
    of arguments. When you are done, you must call shutdown() to stop \
    the running processes. Alternatively, you can use the MultiProcessor \
    as a context manager using a with statement.
    """

    def __init__(self, processor_class: Type, proc_init_arguments: Any = None, num_processes: int = None) -> None:
        """
        Create a MultiProcessor.

        Args:
            processor_class: Class to be instantiated in each process.
            proc_init_arguments: Arguments to pass to the constructor. \
                    Either a single object or a tuple of multiple \
                    arguments, which will be unpacked.
            num_processes: Number of processes to use. Defaults to the \
                    number of CPUs in the machine.
        """
        if not proc_init_arguments:
            proc_init_arguments = ()

        if not isinstance(proc_init_arguments, tuple):
            proc_init_arguments = (proc_init_arguments,)

        if not num_processes:
            num_processes = mp.cpu_count()

        queue_length = 2 * num_processes
        self.__in_queue = mp.Queue(queue_length)
        self.__out_queue = mp.Queue(queue_length)

        self.__processes = []
        for i in range(num_processes):
            process = mp.Process(
                    target=MultiProcessor.__runner,
                    args=(self.__in_queue, self.__out_queue, processor_class, proc_init_arguments))
            process.start()
            self.__processes.append(process)


    def __enter__(self) -> 'MultiProcessor':
        return self


    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.shutdown()


    def processed(self, items: List, progress_label: str=None) -> Generator:
        """
        Processes the given items in random order.

        Sends each item in items to a processor, yielding results in \
        random order. Processing is done in parallel. If an item is a \
        tuple, then it will be unpacked and sent to the processor as \
        multiple arguments for its process() function.

        Args:
            items: A list of items to process.
            progress_label: Label for printed progress information. \
                    Set to None or omit to not print progress.

        Returns:
            A list of results, in arbitrary order.
        """
        num_queued = 0
        num_results = 0
        while num_results < len(items):
            while num_queued < len(items) and not self.__in_queue.full():
                next_item = items[num_queued]
                if not isinstance(next_item, tuple):
                    next_item = (next_item,)

                self.__in_queue.put(next_item)
                num_queued += 1

            if num_results < len(items):
                yield self.__out_queue.get()
                num_results += 1
                if progress_label is not None:
                    print_progress(num_results, len(items), progress_label)

    def shutdown(self) -> None:
        """
        Shut down the MultiProcessor.

        Calling this function stops the external processes, and frees \
        the associated resources.
        """
        for _ in self.__processes:
            self.__in_queue.put(_Stop())

        for proc in self.__processes:
            proc.join()


    @staticmethod
    def __runner(in_queue: mp.Queue, out_queue: mp.Queue, processor_class: Type, proc_init_arguments: Tuple) -> None:
        processor = processor_class(*proc_init_arguments)

        while True:
            arguments = in_queue.get()
            if isinstance(arguments, _Stop):
                break
            result = processor.process(*arguments)
            out_queue.put(result)

        sys.exit()
