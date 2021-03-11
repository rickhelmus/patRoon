function getAnnotationPageElement()
{
    return document.getElementById("feature-groups").parentElement.parentElement;
}

function showPlot(group, cl, paths)
{
    el = document.getElementById(cl);
    el.style.display = "flex";
    el.style.maxWidth = "100%";
    el.src = paths[parseInt(group) - 1];
}

function showAnnotation(group, type)
{
    const annElements = document.getElementsByClassName('annotationClass');
    for (var i=0; i<annElements.length; i++)
        annElements[i].style.display = (annElements[i].classList.contains(type)) ? 'flex' : 'none';

    showPlot(group, "EICAnn", chromPaths);

    qu = "#" + type + "Table .dataTable";
    $(qu).DataTable().column(0).search("^" + group + "$", true, false).draw();
    $(qu).DataTable().columns.adjust().draw();

    selectFGroupRow(group, "fGroupsTable");
    
    document.getElementById("noAnnotationSelected").style.display = 'none';
}

function showCompoundsCluster(group)
{
    const type = "compscl_ann-" + group;
    const annElements = document.getElementsByClassName('annotationClass');
    for (var i=0; i<annElements.length; i++)
        annElements[i].style.display = (annElements[i].classList.contains(type)) ? 'flex' : 'none';

    showPlot(group, "EICAnn", chromPaths);
    selectFGroupRow(group, "fGroupsTable");
    document.getElementById("noAnnotationSelected").style.display = 'none';
}

function disableAllAnnotations(cl)
{
    var comps = document.getElementsByClassName(cl);
    for (var i=0; i<comps.length; i++)
        comps[i].style.display = 'none';
    //$(".dataTable").DataTable().columns.adjust().draw(); // fixup feature group table
}

function selectFGroupRow(group, cl)
{
    var table = $("#" + cl + " .dataTable").DataTable();
    
    table.$('tr.selected').removeClass('selected'); // remove any current selections
    
    var indexes = table.rows().eq(0).filter(function(ind)
    {
        return table.cell(ind, 0).data() == group ? true : false;
    });
    table.rows(indexes).nodes().to$().addClass("selected");
}

function initAnnotation()
{
    disableAllAnnotations("annotationClass");
    $('.dataTable').DataTable().columns.adjust().draw();
}

function initTPs()
{
    disableAllAnnotations("TPsClass");
    $('.dataTable').DataTable().columns.adjust().draw();
}

function showTPs(cmp, group)
{
    const elements = document.getElementsByClassName("TPsClass");
    for (var i=0; i<elements.length; i++)
        elements[i].style.display = 'flex';

    showPlot(group, "precEIC", chromPaths);
    showPlot(cmp, "precStruct", precStructPaths); // UNDONE: this may not be there
    showPlot(cmp, "precInt", precIntPaths);

    qu = "#TPsTable .dataTable";
    $(qu).DataTable().column(0).search("^" + group + "$", true, false).draw();
    $(qu).DataTable().columns.adjust().draw();

    selectFGroupRow(group, "precursorsTable");
    document.getElementById("noPrecSelected").style.display = 'none';
}

$(document).ready(function()
{
    //setTimeout(disableAllAnnotations, 5000); // HACK: wait a bit so that HTML instances are available
});
