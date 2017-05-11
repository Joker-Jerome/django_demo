from django import forms
from .forms import DocumentForm
from django.shortcuts import render
from django.shortcuts import redirect
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from .centrality import calculate

def model_form_upload(request):
    if request.method == 'POST':
        form = DocumentForm(request.POST, request.FILES)
        if form.is_valid():
            form.save()
            file_name = "media/datasets/" + str(request.FILES["document"])
            calculate(file_name)
            #return render(request,'upload/csv-to-html-table/index.html')
            #return redirect("blog")
            return redirect("result")
    else:
        form = DocumentForm()
    return render(request, 'upload/model_form_upload.html', {
                  'form': form
                  })
def showresult(request):
    return render(request,'upload/result.html')

def showresult_snp(request):
    return render(request,'upload/csv-to-html-table/index_snp.html')

def showresult_nonsnp(request):
    return render(request,'upload/csv-to-html-table/index_nonsnp.html')


