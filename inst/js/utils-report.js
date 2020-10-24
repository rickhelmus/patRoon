function getAnnotationPageElement()
{
    return document.getElementById("feature-groups").parentElement.parentElement;
}

function showEICs(group)
{
    el = document.getElementById("EIC");
    el.style.display = "flex";
    el.src = EICPaths[parseInt(group) - 1];
    document.getElementById("noAnnotationSelected").style.display = 'none';
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
    var table = $("#fGroupsTable .dataTable").DataTable();
    
    table.$('tr.selected').removeClass('selected'); // remove any current selections
    
    var indexes = table.rows().eq(0).filter(function(ind)
    {
        return table.cell(ind, 0).data() == group ? true : false;
    });
    table.rows(indexes).nodes().to$().addClass("selected");
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
