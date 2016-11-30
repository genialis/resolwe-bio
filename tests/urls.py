from __future__ import absolute_import, division, print_function, unicode_literals

from django.conf.urls import include, url

from rest_framework import routers

from resolwe_bio.views import SampleViewSet, PresampleViewSet
from resolwe_bio.kb.views import FeatureSearchViewSet, FeatureAutocompleteViewSet, FeatureViewSet, MappingViewSet

api_router = routers.DefaultRouter(trailing_slash=False)  # pylint: disable=invalid-name
api_router.register(r'sample', SampleViewSet)
api_router.register(r'presample', PresampleViewSet, 'presample')
api_router.register(r'kb/feature/admin', FeatureViewSet)
api_router.register(r'kb/feature/search', FeatureSearchViewSet, 'kb_feature_search')
api_router.register(r'kb/feature/autocomplete', FeatureAutocompleteViewSet, 'kb_feature_autocomplete')
api_router.register(r'kb/mapping', MappingViewSet)

urlpatterns = [  # pylint: disable=invalid-name
    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),
    url(r'^', include(api_router.urls, namespace='resolwebio-api')),
]
