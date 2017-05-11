from django import forms
from .forms import DocumentForm
from django.shortcuts import render
from django.shortcuts import redirect
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from .deleteriousness_estimation_v2 import calculate

def model_form_upload(request):
    if request.method == 'POST':
        form = DocumentForm(request.POST, request.FILES)
        if form.is_valid():
            form.save()
            file_name = "media/datasets/" + str(request.FILES["document"])
            calculate(file_name)
            #return render(request,'deleteriousness_estimation/csv-to-html-table/index.html')
            #return redirect("blog")
            return redirect("ds")
    else:
        form = DocumentForm()
    return render(request, 'deleteriousness_estimation/model_form_upload.html', {
                  'form': form
                  })

def showresult(request):
    return render(request,'deleteriousness_estimation/csv-to-html-table/index.html')
