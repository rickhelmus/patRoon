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
    
    qu = "#" + type + "Table .dataTable";
    $(qu).DataTable().column(0).search(group).draw();
    $(qu).DataTable().columns.adjust().draw();
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
    $(".dataTable").DataTable().columns.adjust().draw(); // fixup feature group table
}

$(document).ready(function()
{
    //setTimeout(disableAllAnnotations, 5000); // HACK: wait a bit so that HTML instances are available
});
