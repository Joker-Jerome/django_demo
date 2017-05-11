from django.conf.urls import url,include
from . import views

urlpatterns = [
               url(r'^$', views.model_form_upload,name='upload'),
               url(r'^result$',views.showresult,name='result'),
               url(r'^snp_results/$',views.showresult_snp,name='snp_result'),
               url(r'^nonsnp_results/$',views.showresult_nonsnp,name='nonsnp_result'),

               
               ]

