function getFeatTableIDs()
{
    const opts = document.getElementById("view-select").options;
    return Array.from(opts).map(o => "detailsTab" + o.value);
}

function getFeatTableElements()
{
    return getFeatTableIDs().map(id => document.getElementById(id));
}

function getFeatTableInstances()
{
    return getFeatTableIDs().map(id => Reactable.getInstance(id));
}

function updateView(sel)
{
    sel = "detailsTab" + sel;
    getFeatTableElements().forEach(el => el.style.display = (el.id === sel) ? "" : "none");
}

function showFeatCols(column, show)
{
    const tabIDs = getFeatTableIDs();
    debugger;
    if (column === "chrom_large")
        tabIDs.forEach(id => Reactable.toggleHideColumn(id, column, show));
}

function showTPGraph(cmp)
{
    TPGraphs = document.querySelectorAll('[id ^= "TPGraph_"]');
    const elName = "TPGraph_" + cmp
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
}
