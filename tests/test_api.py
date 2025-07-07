import unittest
from unittest.mock import patch, MagicMock
from phytominer import api

class TestApiModule(unittest.TestCase):

    @patch('requests.get')
    def test_fetch_gene_data_success(self, mock_get):
        # Arrange
        expected_response = {'gene_id': 'AT1G01110', 'description': 'Sample Gene'}
        mock_get.return_value = MagicMock(status_code=200)
        mock_get.return_value.json.return_value = expected_response

        result = api.fetch_gene_data('AT1G01110')

        # Assert
        self.assertEqual(result, expected_response)
        mock_get.assert_called_once_with('https://api.example.com/genes/AT1G01110')

    @patch('requests.get')
    def test_fetch_gene_data_failure(self, mock_get):
        # Arrange
        mock_get.return_value = MagicMock(status_code=404)
        mock_get.return_value.json.return_value = {'error': 'Not found'}
        # Act
        result = api.fetch_gene_data('INVALID_ID')
        # Assert
        self.assertIsNone(result)
        mock_get.assert_called_once_with('https://api.example.com/genes/INVALID_ID')

if __name__ == '__main__':
    unittest.main()
