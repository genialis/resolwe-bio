# pylint: disable=missing-docstring
from unittest import TestCase

import six


class TestTools(TestCase):

    # pylint: disable=import-error,no-name-in-module
    def test_tools_is_not_a_package(self):
        if six.PY2:
            with self.assertRaises(ImportError):
                import resolwe_bio.tools
        else:
            # Python 3.3+ introduced implicit namespace packages in PEP 420,
            # therefore importing will not fail, however, namespace packages
            # have no '__file__' attribute
            import resolwe_bio.tools
            with self.assertRaises(AttributeError):
                resolwe_bio.tools.__file__  # pylint: disable=no-member,pointless-statement
