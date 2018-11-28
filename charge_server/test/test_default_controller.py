# coding: utf-8

from __future__ import absolute_import

from flask import json
from six import BytesIO

from charge_server.test import BaseTestCase


class TestDefaultController(BaseTestCase):
    """DefaultController integration test stubs"""

    def test_charge_molecule(self):
        """Test case for charge_molecule

        Submit a molecule for charging
        """
        body = 'body_example'
        response = self.client.open(
            '/charge_assign',
            method='POST',
            data=json.dumps(body),
            content_type='text/plain')
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))


if __name__ == '__main__':
    import unittest
    unittest.main()
