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
    else if (which === "ISTDs")
        el = document.getElementById("struct_view-istd");
    else if (which === "Components")
        el = document.getElementById("chrom_view-component");
    else if (which === "TPs")
        el = document.getElementById("int_plot-parent");

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
    
    if (document.getElementById("suspAnnTab"))
    {
        const showSusps = sel === "Suspects" || sel === "TPs";
        
        showFeatureTab("Suspect annotation", showSusps);
        
        for (ann of [ "formulas", "compounds" ])
        {
            if (document.getElementById(ann + "Tab"))
            {
                let suspCheckEl = document.getElementById(ann + "-susp_only");
                const d = (showSusps) ? "" : "none"
                suspCheckEl.style.display = d;
                Array.from(suspCheckEl.labels).forEach(l => l.style.display = d);
                if (!showSusps)
                    toggleAnnOnlySusp(ann, false); // UNDONE: restore selection when going back to suspect view?
            }
        }
    }
    if (getViews().includes("TPs") && document.getElementById("similarity_spec"))
        showFeatureTab("Parent similarity", sel === "TPs");
    
    const r = Reactable.getInstance(tid).rowsById[Reactable.getState(tid).meta.selectedRow];
    updateFeatTabRowSel(r.values, r.index);
}

function updateFeatTabRowSel(rowValues, rowIndex)
{
    if (rowIndex == undefined)
        return; // don't process clicks of parent rows 
    
    const tabEl = getSelFGTableElement();
    const grp = rowValues.group;
    const plots = Reactable.getState(tabEl).meta.plots;
    
    Reactable.setMeta(tabEl, { selectedRow: rowIndex });

    Reactable.setFilter('featuresTab', 'group', grp);
    document.getElementById('int_plot').src = plots.intPlots[grp];
    if (document.getElementById('spectrumMS'))
    {
        document.getElementById('spectrumMS').src = plots.MSPeakLists[grp].MS;
        Reactable.setFilter('MSPLTab', 'group', grp);
        document.getElementById('spectrumMSMS').src = plots.MSPeakLists[grp].MSMS;
        Reactable.setFilter('MSMSPLTab', 'group', grp);
    }
    if (document.getElementById('formulasTab'))
        Reactable.setFilter('formulasTab', 'group', grp);
    if (document.getElementById('compoundsTab'))
    {
        Reactable.setFilter('compoundsTab', 'group', grp);
        document.getElementById('openMF').href = Reactable.getState('compoundsTab').meta.mfWebLinks[grp];
    }
    
    const ccd = document.getElementById('comps_cluster-dendro');
    if (ccd)
    {
        ccd.src = plots.compsCluster[grp].dendro;
        Array.from(document.getElementsByClassName('mcs')).forEach(el => el.style.display = (el.classList.contains('mcs-' + grp)) ? '' : 'none');
    }

    if (tabEl === "detailsTabSuspects")
    {
        const structEl = document.getElementById('struct_view-suspect');
        structEl.src = plots.structs[rowValues.susp_InChIKey] || "";
        Reactable.setFilter('suspInfoTab', 'name', rowValues.susp_name);
        if (document.getElementById('suspAnnTab'))
            Reactable.setFilter('suspAnnTab', 'suspID', rowValues.susp_name + '-' + rowValues.group);
    }
    else if (tabEl === "detailsTabISTDs")
    {
        const structEl = document.getElementById('struct_view-istd');
        structEl.src = plots.structs[rowValues.InChIKey] || "";
        Reactable.setFilter('ISTDInfoTab', 'name', rowValues.susp_name);
    }
    else if (tabEl === "detailsTabComponents")
    {
        let chromEl = document.getElementById('chrom_view-component');
        let specEl = document.getElementById('spectrum_view-component');
        let profileRelEl = document.getElementById('profileRel_view-component');
        let profileAbsEl = document.getElementById('profileAbs_view-component');
        const pl = plots.components.components[rowValues.component];
        chromEl.src = pl.chrom;
        specEl.src = pl.spec;
        if (profileRelEl != undefined)
        {
            profileRelEl.src = pl.profileRel;
            profileAbsEl.src = pl.profileAbs;
        }
        if (document.getElementById('componentInfoTab'))
            Reactable.setFilter('componentInfoTab', 'name', rowValues.component);
    }
    else if (tabEl === "detailsTabTPs")
    {
        const chromEl = document.getElementById('chrom_view-tp');
        const intEl = document.getElementById('int_plot-parent');
        const specSimEl = document.getElementById('similarity_spec');
        
        if (chromEl)
            chromEl.src = plots.chromsLarge[rowValues.parent_group];
        intEl.src = plots.intPlots[rowValues.parent_group];
        
        if (document.getElementById('parentInfoTab'))
        {
            document.getElementById('struct_view-parent').src = plots.structs[rowValues.parent_susp_InChIKey];
            Reactable.setFilter('parentInfoTab', 'name', rowValues.parent_susp_name);
            document.getElementById('struct_view-tp').src = plots.structs[rowValues.susp_InChIKey];
            Reactable.setFilter('TPInfoTab', 'name', rowValues.susp_name);
        }
        
        if (plots.TPs.length > 0)
        {
            specSimEl.src = plots.TPs[rowValues.component][rowValues.cmpIndex];
            specSimEl.style.display = ''; // may have been hidden if a previous img didn't exist
        }
        if (document.getElementById('suspAnnTab'))
            Reactable.setFilter('suspAnnTab', 'suspID', rowValues.susp_name + '-' + rowValues.group);
        Reactable.setFilter('similarityTab', 'cmpID', rowValues.component + '-' + rowValues.cmpIndex);
        
        showTPGraph(rowValues.component);
    }
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

function showFeatQualityCols(show)
{
    const cols = Reactable.getState("featuresTab").meta.featQualCols;
    cols.forEach(col => Reactable.toggleHideColumn("featuresTab", col, !show));
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

function filtRangeModalInit(tab, col)
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

function filtSelModalInit(tab, col)
{
    const mname = "filter_" + col;
    const curF = Reactable.getState(tab).meta[mname];

    let lg = document.getElementById("filterSelectListGroup");
    while (lg.firstChild)
        lg.removeChild(lg.firstChild);
    
    // UNDONE: this always adds all (unique) values, including those that are e.g. filtered out due to group selection    
    const vals = new Set(Reactable.getState(tab).data.map(row => row[col]));
    let checkEls = [ ];
    vals.forEach(function(label)
    {
        let inpEl = document.createElement("input");
        inpEl.classList.add("form-check-input", "me-1");
        inpEl.type = "checkbox";
        inpEl.checked = curF == undefined || curF.has(label);
        
        let labEl = document.createElement("label");
        labEl.classList.add("list-group-item");
        labEl.appendChild(inpEl);
        labEl.appendChild(document.createTextNode(label));
        
        lg.appendChild(labEl);
        checkEls.push(inpEl);
    });
    
    document.getElementById("filtSelApply").addEventListener("click", function(e)
    {
        const checked = new Set();
        checkEls.forEach(function(cel)
        {
            if (cel.checked)
                checked.add(cel.nextSibling.textContent);
        });
        Reactable.setFilter(tab, col, checked);
        Reactable.setMeta(tab, { [mname]: checked });
    }, { once: true });
    
    //window.addEventListener('keydown', filtModalKeyHandler);
}

function filtSelModalToggleAll(e)
{
    Array.from(document.getElementById("filterSelectListGroup").children).forEach(el => el.children[0].checked = e);
}

function applyFilterToggle(tab, nonFiltCol, keepGroup = false)
{
    // HACK: setting a filter will toggle the visibility if all filters were enabled/disabled
    // NOTE: use a non filterable column so nothing gets reset
    // NOTE: this somehow will reset all filter values, so restore group filter if needed...
    const grp = (keepGroup) ? Reactable.getInstance("featuresTab").allColumns.find(c => c.id === "group").filterValue : undefined;
    Reactable.setFilter(tab, nonFiltCol, undefined);
    if (keepGroup)
        Reactable.setFilter(tab, "group", grp); 
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
    
    // NOTE: we only have to do this on the active table as the rest will be re-drawn when activated
    applyFilterToggle(getSelFGTableElement(), "chrom_small");
}

function toggleFeatFilters(e)
{
    // as above, for feature table
    Reactable.getInstance("featuresTab").allColumns.forEach(function(col)
    {
        if (col.name !== "chromatogram" && col.name !== "group") // UNDONE: do this more elegantly?
            col.filterable = e;
    })
    
    applyFilterToggle("featuresTab", "group", true);
}

function toggleFormFilters(e)
{
    // as above, for formulas table
    const skipCols = [ ".details", "group", "spectrum", "scorings" ];
    Reactable.getInstance("formulasTab").allColumns.forEach(function(col)
    {
        if (!skipCols.includes(col.id)) // UNDONE: do this more elegantly?
            col.filterable = e;
    })
    applyFilterToggle("formulasTab", "spectrum", true);
}

function toggleCompFilters(e)
{
    // as above, for compounds table
    const skipCols = [ ".details", "group", "structure", "spectrum", "scorings" ];
    Reactable.getInstance("compoundsTab").allColumns.forEach(function(col)
    {
        if (!skipCols.includes(col.id)) // UNDONE: do this more elegantly?
            col.filterable = e;
    })
    applyFilterToggle("compoundsTab", "spectrum", true);
}

function toggleAnnOnlySusp(wh, e)
{
    const tid = (wh === "formulas") ? "formulasTab" : "compoundsTab";
    if (!e)
        Reactable.setFilter(tid, "suspect", undefined);
    else
    {
        const fgTab = getSelFGTableElement();
        const curRow = Reactable.getState(fgTab).meta.selectedRow;
        Reactable.setFilter(tid, "suspect", Reactable.getInstance(fgTab).rowsById[curRow].values.susp_name);
    }
}

function showFeatureTab(tabName, enable)
{
    const featTabsEl = document.getElementById("fGroupSelTabs");
    let tabEl = featTabsEl.querySelectorAll('[data-value="' + tabName + '"]')[0];
    let tabLiEl = tabEl.parentElement;
    if (enable)
    {
        $(tabLiEl).show();
    }
    else
    {
        $(tabLiEl).hide();
        // UNDONE: seems active class is not always applied on the same elements?
        if (tabLiEl.classList.contains("active") || tabEl.classList.contains("active"))
        {
            // Activate features tab, assuming it is always available
            let featTab = featTabsEl.querySelectorAll('[data-value="Features"]')[0];
            $(featTab).tab("show");
            
            // remove active classes, seems to be that this needs to be done last
            tabLiEl.classList.remove("active");
            tabEl.classList.remove("active");
        }
    }
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
