from django.conf.urls import url,include
from . import views

urlpatterns = [
               url(r'^$', views.model_form_upload,name='upload_ds'),
               url(r'^del$',views.showresult,name='ds'),
               
               ]

