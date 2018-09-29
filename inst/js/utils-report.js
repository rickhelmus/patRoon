function getAnnotationPageElement()
{
    return document.getElementById("feature-groups").parentElement.parentElement;
}

function showEICs(group)
{
    const annElements = document.getElementsByClassName('annotationClass');
    const EICId = "EIC-" + group;
    for (var i=0; i<annElements.length; i++)
    {
        if (annElements[i].classList.contains(EICId))
            annElements[i].style.display = 'flex';
    }
    
    document.getElementsByClassName("noAnnotationSelected")[0].style.display = 'none';    
}

function showAnnotation(group, type)
{
    const annElements = document.getElementsByClassName('annotationClass');
    for (var i=0; i<annElements.length; i++)
        annElements[i].style.display = (annElements[i].classList.contains(type)) ? 'flex' : 'none';

    showEICs(group);

    qu = "#" + type + "Table .dataTable";
    $(qu).DataTable().column(0).search("^" + group + "$", true, false).draw();
    $(qu).DataTable().columns.adjust().draw();
    
    selectFGroupRow(group);
}

function showCompoundsCluster(group)
{
    const type = "compscl_ann-" + group;
    const annElements = document.getElementsByClassName('annotationClass');
    for (var i=0; i<annElements.length; i++)
        annElements[i].style.display = (annElements[i].classList.contains(type)) ? 'flex' : 'none';
    
    showEICs(group);
    selectFGroupRow(group);
}

function disableAllAnnotations()
{
    var comps = document.getElementsByClassName('annotationClass');
    for (var i=0; i<comps.length; i++)
        comps[i].style.display = 'none';
    //$(".dataTable").DataTable().columns.adjust().draw(); // fixup feature group table
}

function selectFGroupRow(group)
{
    qu = "#fGroupsTable .dataTable";
    $(qu).DataTable().$('tr.selected').removeClass('selected'); // remove any current selections
    $(qu).DataTable().row(group-1).node().classList.add("selected");
}

function initAnnotation()
{
    disableAllAnnotations();
    $('.dataTable').DataTable().columns.adjust().draw();
}

$(document).ready(function()
{
    //setTimeout(disableAllAnnotations, 5000); // HACK: wait a bit so that HTML instances are available
});
