# coding: utf-8

from __future__ import absolute_import

import os

from flask import json
from six import BytesIO

from charge_server.test import BaseTestCase
from charge_server import charge_server


class TestDefaultController(BaseTestCase):
    """DefaultController integration test stubs"""

    def test_charge_molecule(self):
        """Test case for charge_molecule

        Submit a molecule for charging
        """
        repo = os.path.join(os.path.dirname(__file__), 'test_repo.zip')
        charge_server.init(repo)

        test_lgf = ('@nodes\n'
                    'label\tlabel2\tatomType\tinitColor\t\n'
                    '1\tC1\t12\t0\t\n'
                    '2\tH1\t20\t0\t\n'
                    '3\tH2\t20\t0\t\n'
                    '4\tH3\t20\t0\t\n'
                    '5\tH4\t20\t0\t\n'
                    '@edges\n'
                    '\t\tlabel\t\n'
                    '1\t2\t0\t\n'
                    '1\t3\t1\t\n'
                    '1\t4\t2\t\n'
                    '1\t5\t3\t\n')
        query_string = [('total_charge', 0)]
        response = self.client.open(
            '/charge_assign',
            method='POST',
            data=test_lgf,
            content_type='text/plain',
            query_string=query_string)
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))


if __name__ == '__main__':
    import unittest
    unittest.main()
