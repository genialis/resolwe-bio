"""Routing configuration for Django Channels."""
# from channels.routing import route

from channels.routing import route_class

from resolwe.flow.managers import manager


manager.update_routing()
channel_routing = [  # pylint: disable=invalid-name
    route_class(type(manager)),
]
