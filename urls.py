from django.conf.urls.defaults import *

from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    # Example:
    # (r'^dnafilter_d2/', include('dnafilter_d2.foo.urls')),
    
    #(r'',)
    
    url(r'^$', 'dnafilter.views.index'),
    url(r'^site_media/(?P<path>.*)$', 'django.views.static.serve',
    {'document_root': '/home/sebastian/Projects/dnafilter/dnafilter_d2/dnafilter/static'}),
    url(r'^filter$', 'dnafilter.views.filter'),

    (r'^admin/', include(admin.site.urls)),
)
