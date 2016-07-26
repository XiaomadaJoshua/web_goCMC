from django.conf.urls import url

from . import views

app_name = 'doseCalc'
urlpatterns = [
	url(r'^$', views.IndexView.as_view(), name='index'),
	url(r'^submit/$', views.submit, name='submit'),
	url(r'^calculating/$', views.calculating, name='calculating'),
]
