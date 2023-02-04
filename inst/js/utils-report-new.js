function getViews()
{
    return Array.from(document.getElementById("view-select").options).map(o => o.value);
}

function getFGTableIDs()
{
    return getViews().map(v => "detailsTab" + v);
}

function getFGTableElements()
{
    return getFGTableIDs().map(id => document.getElementById(id));
}

function getFGTableInstances()
{
    return getFGTableIDs().map(id => Reactable.getInstance(id));
}

function getSelFGTableElement()
{
    return "detailsTab" + document.getElementById("view-select").value;
}

function getNavTab(which)
{
    // HACK: the following functions work their way up from a child element to the relevant parent, since it seems no
    // class or ID can be set for bslib::nav_tab_card()

    let el;
    if (which === "Suspects")
        el = document.getElementById("struct_view-suspect");
    else if (which === "Components")
        el = document.getElementById("chrom_view-component");
    else if (which === "TPs")
        el = document.getElementById("chrom_view-tp");

    for (var i=0; i<6; i++)
        el = el.parentElement;
    
    return el;
}

function updateView(sel)
{
    tid = "detailsTab" + sel;
    getFGTableElements().forEach(el => el.style.display = (el.id === tid) ? "" : "none");
    document.getElementById("fg-expand").style.display = (sel !== "Plain") ? "" : "none";
    
    document.getElementsByClassName("bottomLayout")[0].style["grid-template-columns"] = (sel === "Plain") ? "1fr" : "1fr 2fr";
    
    getViews().forEach(function(v)
    {
        if (v === "Plain")
            return;
        let el = getNavTab(v);
        el.classList.toggle("d-none", v !== sel)
    })
}

function showFGCols(column, show)
{
    const tabIDs = getFGTableIDs();
    tabIDs.forEach(function(id)
    {
        const cols = Reactable.getState(id).meta.colToggles[column];
        if (Array.isArray(cols))
            cols.forEach(col => Reactable.toggleHideColumn(id, col, !show));
        else
            Reactable.toggleHideColumn(id, cols, !show);
    })
    if (column === "chrom_large")
        tabIDs.forEach(id => Reactable.toggleHideColumn(id, "chrom_small", show));
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

function initRangeModal(tab, col)
{
    const mname = "filter_" + col;
    const curF = Reactable.getState(tab).meta[mname];
    
    document.getElementById("filtNumMin").value = (curF == undefined) ? "" : curF[0];
    document.getElementById("filtNumMax").value = (curF == undefined) ? "" : curF[1];
    
    document.getElementById("filtNumApply").addEventListener("click", function(e)
    {
        const r = [ document.getElementById("filtNumMin").value, document.getElementById("filtNumMax").value ];
        Reactable.setFilter(tab, col, r);
        Reactable.setMeta(tab, { [mname]: r });
    }, { once: true });
    
    window.addEventListener('keydown', filtModalKeyHandler);
}

function filtModalKeyHandler(e)
{
    // based on https://stackoverflow.com/a/41055853
    if ($("#filterRangeModal").hasClass("show") && (e.keycode == 13 || e.which == 13))
    {
        document.getElementById("filtNumApply").click();
        window.removeEventListener('keydown', filtModalKeyHandler);    
    }
}

function filtNumClear()
{
    document.getElementById("filtNumMin").value = "";
    document.getElementById("filtNumMax").value = "";
    document.getElementById("filtNumApply").click();
}

function toggleFGFilters(e)
{
    const tabIDs = getFGTableIDs();
    tabIDs.forEach(function(id)
    {
        Reactable.getInstance(id).allColumns.forEach(function(col)
        {
            if (col.name !== "chromatogram") // UNDONE: do this more elegantly?
                col.filterable = e;
        })
    })
    // HACK: this also redraws table to toggle the filter row
    // NOTE: we only have to do this on the active table as the rest will be re-drawn when activated
    // NOTE: use a non filterable column so nothing gets reset
    Reactable.setFilter(getSelFGTableElement(), "chrom_small", undefined);
}

$(document).ready(function() {
    // Image zooming, based on https://stackoverflow.com/a/57694495
    $('body').prepend("<div class=\"zoomDiv\"><img src=\"\" class=\"zoomImg\"></div>");
    $('body').on('click', 'img:not(.zoomImg, .noZoomImg)', function() {
       $('.zoomImg').attr('src', $(this).attr('src'));
       $('.zoomDiv').css({opacity: '1', width: '70%'});
    });
    $('img.zoomImg').click(function() {
       $('.zoomDiv').css({opacity: '0', width: '0%'});
    });
    
    updateView("Plain");
    
    document.getElementById("filterRangeModal").addEventListener('shown.bs.modal', () => {
        document.getElementById("filtNumMin").focus();
    });
});
