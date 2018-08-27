function getAnnotationPageElement()
{
    return document.getElementById("feature-groups").parentElement.parentElement;
}

function showAnnotation(group, type)
{
    clearComponentSpecs();
    
    const annElements = document.getElementsByClassName('annotationClass');
    const EICId = "EIC-" + group;
    for (var i=0; i<annElements.length; i++)
    {
        annElements[i].style.display = (annElements[i].classList.contains(type) ||
                                        annElements[i].classList.contains(EICId)) ? 'flex' : 'none';
    }
    
    document.getElementsByClassName("noAnnotationSelected")[0].style.display = 'none';
    
    if (type == "compounds")
        $('#compoundsTable .dataTable').DataTable().column(0).search(group).draw();
    else
        $('#formulasTable .dataTable').DataTable().column(0).search(group).draw();
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
    //setTimeout(function() { $('.annotationClass').hide(); }, 1500);
    setTimeout(disableAllAnnotations, 5); // HACK: wait a bit so that HTML instances are available
});
