from django.urls import include, path
from rest_framework import routers

from resolwe.api_urls import api_router as resolwe_router
from resolwe.flow.views import EntityViewSet

from resolwe_bio.filters import BioEntityFilter
from resolwe_bio.kb.views import FeatureViewSet, MappingSearchViewSet
from resolwe_bio.variants.views import (
    VariantAnnotationViewSet,
    VariantCallViewSet,
    VariantExperimentViewSet,
    VariantViewSet,
)

from .routers import SearchRouter

EntityViewSet.filterset_class = BioEntityFilter

api_router = routers.DefaultRouter(trailing_slash=False)
api_router.register(r"sample", EntityViewSet)
api_router.register(r"variant", VariantViewSet)
api_router.register(r"variant_annotations", VariantAnnotationViewSet)
api_router.register(r"variant_calls", VariantCallViewSet)
api_router.register(r"variant_experiment", VariantExperimentViewSet)

search_router = SearchRouter(trailing_slash=False)
search_router.register(r"kb/feature", FeatureViewSet, "kb_feature")
search_router.register(r"kb/mapping/search", MappingSearchViewSet, "kb_mapping_search")

urlpatterns = [
    path("api-auth/", include("rest_framework.urls", namespace="rest_framework")),
    # XXX: Temporary fix to work with Resolwe 2.0.0, which requires 'resolwe-api' namespace to be available when
    # reporting errors when running processes.
    path("api-resolwe/", include((resolwe_router.urls, "resolwe-api"))),
    path(
        "api/",
        include(
            (
                api_router.urls + search_router.urls + resolwe_router.urls,
                "resolwebio-api",
            )
        ),
    ),
]
