// SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
//
// SPDX-License-Identifier: GPL-3.0-only

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

    filterDTRows(group, type + "Table");
    selectDTRow(group, "fGroupsTable");
    
    document.getElementById("noAnnotationSelected").style.display = 'none';
}

function showCompoundsCluster(group)
{
    const type = "compscl_ann-" + group;
    const annElements = document.getElementsByClassName('annotationClass');
    for (var i=0; i<annElements.length; i++)
        annElements[i].style.display = (annElements[i].classList.contains(type)) ? 'flex' : 'none';

    showPlot(group, "EICAnn", chromPaths);
    selectDTRow(group, "fGroupsTable");
    document.getElementById("noAnnotationSelected").style.display = 'none';
}

function disableAllAnnotations(cl)
{
    var comps = document.getElementsByClassName(cl);
    for (var i=0; i<comps.length; i++)
        comps[i].style.display = 'none';
    //$(".dataTable").DataTable().columns.adjust().draw(); // fixup feature group table
}

function filterDTRows(what, cl)
{
    qu = "#" + cl + " .dataTable";
    $(qu).DataTable().column(0).search("^" + what + "$", true, false).draw();
    $(qu).DataTable().columns.adjust().draw();
}

function selectDTRow(index, cl, cell = null)
{
    var table = $("#" + cl + " .dataTable").DataTable();
    
    table.$((cell) ? "td.selected" : "tr.selected").removeClass("selected"); // remove any current selections
    
    var indexes = table.rows().eq(0).filter(function(ind)
    {
        return table.cell(ind, 0).data() == index ? true : false;
    });
    
    if (cell)
        table.cells(indexes, cell).nodes().to$().addClass("selected");
    else
        table.rows(indexes).nodes().to$().addClass("selected");
}

function initAnnotation()
{
    disableAllAnnotations("annotationClass");
    $('.dataTable').DataTable().columns.adjust().draw();
}

function initTPs()
{
    filterDTRows("nothing", "parentsPlotsTable");
    filterDTRows("nothing", "TPsTable");
    
    var table = $("#parentsTable .dataTable").DataTable();
    var IDs = Array.from(table.column(0).data());
    for (const i of IDs)
    {
        el = document.getElementById("TPGraph_" + i);
        el.style.display = "none";
    }
}

function showTPs(cmp, group)
{
    filterDTRows(cmp, "parentsPlotsTable");
    filterDTRows(cmp, "TPsTable");
    selectDTRow(cmp, "parentsTable", 1);
    
    TPGraphs = document.querySelectorAll('[id ^= "TPGraph_"]');
    elName = "TPGraph_" + cmp
    for (var i=0; i<TPGraphs.length; i++)
    {
        if (TPGraphs[i].id == elName && TPGraphs[i].children.length > 0) // NOTE: no children if plot couldn't be made
        {
            TPGraphs[i].style.display = "";
            document.getElementById("graphTPGraph_" + cmp).chart.fit(); // HACK: reset zoom
        }
        else
            TPGraphs[i].style.display = "none";
    }
    
    $("#parentsTable .dataTable").DataTable().columns.adjust().draw(false);
}

$(document).ready(function()
{
    //setTimeout(disableAllAnnotations, 5000); // HACK: wait a bit so that HTML instances are available
});
