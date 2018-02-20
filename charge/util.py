
def print_progress(iteration:int,
                   total:int,
                   prefix:str='',
                   suffix:str='',
                   decimals:int=1,
                   length:int=100,
                   fill:str='#') -> None:
    """Call in a loop to create terminal progress bar
    
        :param iteration: current iteration
        :type iteration: int
        :param total: total iterations
        :type total: int
        :param prefix: prefix string
        :type prefix: str
        :param suffix: suffix string
        :type suffix: str
        :param decimals: positive number of decimals in percent complete
        :type decimals: int
        :param length: character length of bar
        :type length: int
        :param fill : bar fill character
        :type fill: str
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '', flush=True)

    if iteration == total:
        print()