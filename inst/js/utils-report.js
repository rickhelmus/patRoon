function getAnnotationPageElement()
{
    return document.getElementById("feature-groups").parentElement.parentElement;
}

function showAnnotation(group, type)
{
    $('#compoundsTable .dataTable').DataTable().column(0).search(group).draw();
    return;
    
    clearComponentSpecs();
    
    var comps = document.getElementsByClassName('annotationClass');
    var annotationId = type + "_" + "ann-" + group;
    var EICId = "EIC-" + group;
    for (var i=0; i<comps.length; i++)
        comps[i].style.display = (comps[i].classList.contains(annotationId) || comps[i].classList.contains(EICId)) ? 'flex' : 'none';

    document.getElementById("noAnnotationSelected").style.display = 'none';
}

function showComponentSpec(component, parentID)
{
    // component specs are treated differently as they are cloned
    clearComponentSpecs();
    disableAllAnnotations();
    
    var comp = document.getElementById(component).parentElement;
    var clonedComp = comp.cloneNode(true);
    clonedComp.id = clonedComp.id + "-clone";
    clonedComp.className += " clonedComp";
    
    getAnnotationPageElement().appendChild(clonedComp);
}

function clearComponentSpecs()
{
    $('.clonedComp').remove();
}

function disableAllAnnotations()
{
    var comps = document.getElementsByClassName('annotationClass');
    for (var i=0; i<comps.length; i++)
        comps[i].style.display = 'none';
}

$(document).ready(function()
{
    //     setTimeout(disableAllSpecs, 500); // HACK: wait a bit so that HTML instances are available
});
