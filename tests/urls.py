from __future__ import absolute_import, division, print_function, unicode_literals

from django.conf.urls import include, url

from rest_framework import routers

from resolwe.api_urls import api_router as resolwe_router
from resolwe.elastic import routers as search_routers
from resolwe.flow.views import EntityViewSet
from resolwe_bio.kb.views import (FeatureSearchViewSet, FeatureAutocompleteViewSet, FeatureViewSet,
                                  MappingViewSet, MappingSearchViewSet)

api_router = routers.DefaultRouter(trailing_slash=False)  # pylint: disable=invalid-name
api_router.register(r'sample', EntityViewSet)
api_router.register(r'kb/feature/admin', FeatureViewSet)
api_router.register(r'kb/mapping/admin', MappingViewSet)

search_router = search_routers.SearchRouter(trailing_slash=False)  # pylint: disable=invalid-name
search_router.register(r'kb/feature/search', FeatureSearchViewSet, 'kb_feature_search')
search_router.register(r'kb/feature/autocomplete', FeatureAutocompleteViewSet, 'kb_feature_autocomplete')
search_router.register(r'kb/mapping/search', MappingSearchViewSet, 'kb_mapping_search')

urlpatterns = [  # pylint: disable=invalid-name
    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),
    # XXX: Temporary fix to work with Resolwe 2.0.0, which requires 'resolwe-api' namespace to be available when
    # reporting errors when running processes.
    url(r'^api-resolwe/', include(resolwe_router.urls, namespace='resolwe-api')),
    url(r'^api/', include(api_router.urls + search_router.urls + resolwe_router.urls, namespace='resolwebio-api')),
]
