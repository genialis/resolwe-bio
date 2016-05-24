from __future__ import absolute_import, division, print_function, unicode_literals

from django.conf.urls import include, url

from rest_framework import routers

from resolwe_bio.views import SampleViewSet, PresampleViewSet


api_router = routers.DefaultRouter(trailing_slash=False)  # pylint: disable=invalid-name
api_router.register(r'sample', SampleViewSet)
api_router.register(r'presample', PresampleViewSet, 'presample')

urlpatterns = [  # pylint: disable=invalid-name
    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),
    url(r'^', include(api_router.urls, namespace='resolwebio-api')),
]
