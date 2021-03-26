from unittest import TestCase


class TestTools(TestCase):
    def test_tools_is_not_a_package(self):
        # Python 3.3+ introduced implicit namespace packages in PEP 420,
        # therefore importing will not fail, however, namespace packages
        # have no '__file__' attribute
        import resolwe_bio.tools

        with self.assertRaises(AttributeError):
            result = resolwe_bio.tools.__file__
            if result is None:
                raise AttributeError("Namespace package")
