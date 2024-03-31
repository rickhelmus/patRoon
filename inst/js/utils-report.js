function getViews()
{
    let ret = Array.from(document.getElementById("view-select").options).map(o => o.value);
    if (ret.includes("TPs"))
        ret.push("TPsParents"); // HACK: add parents sub-view of TP view
    return ret;
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
    let sel = document.getElementById("view-select").value;
    if (sel === "TPs" && document.getElementById("viewTPDetailsParents").checked)
        sel = "TPsParents";
    return "detailsTab" + sel;
}

function getNavTab(which)
{return undefined; // UNDONE
    // HACK: the following functions work their way up from a child element to the relevant parent, since it seems no
    // class or ID can be set for bslib::nav_tab_card()

    let el;
    if (which == "Plain")
        el = undefined;
    else if (which === "Suspects")
        el = document.getElementById("struct_view-suspect");
    else if (which === "ISTDs")
        el = document.getElementById("struct_view-istd");
    else if (which === "Components")
        el = document.getElementById("chrom_view-component");
    else if (which === "TPs")
    {
        // For TPs it's a bit more complicated as all the tabs of the panel are optional
        // NOTE: we don't have to check TP graphs, as struct_view-parent will then always be present
        el = document.getElementById("chrom_view-parent") ||
             document.getElementById("struct_view-parent") ||
             document.getElementById("int_plot-parent");
    }

    if (el != undefined)
    {
        while (!el.parentElement.classList.contains("bottomLayout"))
            el = el.parentElement;
    }
    
    return el;
}

function setDetailsTablesRatio(fr1, fr2)
{
    let el = document.getElementById("detailsTablesLayout");
    el.style["grid-template-rows"] = `${fr1}fr ${fr2}fr`;
}

function updateView(sel)
{
    if (sel === "TPs" && document.getElementById("viewTPDetailsParents").checked)
        sel = "TPsParents"; // HACK: parents sub-view
    
    tid = "detailsTab" + sel;
    getFGTableElements().forEach(el => el.style.display = (el.id === tid) ? "" : "none");
    
    // UNDONE
    //document.getElementById("fg-expand").style.display = (sel !== "Plain") ? "" : "none";
    
    // UNDONE
    // document.getElementsByClassName("bottomLayout")[0].style["grid-template-columns"] = (sel === "Plain") ? "1fr" : "1fr 2fr";
    
    getViews().forEach(function(v)
    {
        let el = getNavTab(v);
        if (el)
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
    
    const isTPsOrParsView = (sel === "TPs" || sel === "TPsParents");
    let tpc = document.getElementById("TPCandidatesCard");
    if (tpc)
    {
        // NOTE: we need to toggle the visibility of wrapped containers for cards!
        tpc.parentElement.classList.toggle("d-none", sel !== "TPs");
        document.getElementById("parentCard").parentElement.classList.toggle("d-none", !isTPsOrParsView);
        
        document.getElementById("detailsTopRowLayout").style["grid-template-columns"] = (sel === "TPs") ? "1fr 1fr" : "1fr";
    }

    document.getElementById("TPsParBtGrp").classList.toggle("d-none", !isTPsOrParsView);    
    document.getElementById("detailsLayout").style["grid-template-columns"] = (isTPsOrParsView) ? "1fr 5fr" : "1fr";
    
    const r = Reactable.getInstance(tid).rowsById[Reactable.getState(tid).meta.selectedRow];
    updateFeatTabRowSel(r.values, r.index);
}

function updateFeatTabRowSel(rowValues, rowIndex)
{
    if (rowIndex == undefined)
        return; // don't process clicks of parent rows 
    
    const tabEl = getSelFGTableElement();
    const grp = rowValues.group;
    
    Reactable.setMeta(tabEl, { selectedRow: rowIndex });

    Reactable.setFilter('featuresTab', 'group', grp);
    
    if (document.getElementById('concsTab'))
        Reactable.setFilter('concsTab', 'group', grp);
    if (document.getElementById('toxTab'))
        Reactable.setFilter('toxTab', 'group', grp);
    
    let intEl = document.getElementById('int_plot');
    if (intEl)
        intEl.src = reportPlots.intPlots[grp];
    
    if (document.getElementById('MSPLTab'))
    {
        Reactable.setFilter('MSPLTab', 'group', grp);
        Reactable.setFilter('MSMSPLTab', 'group', grp);
        let specEl = document.getElementById('spectrumMS');
        if (specEl) // not present if !settings$MSPeakLists$spectra
        {
            // NOTE: undefined if filtered away
            specEl.src = (reportPlots.MSPeakLists[grp] && reportPlots.MSPeakLists[grp].MS) || "";
            document.getElementById('spectrumMSMS').src = (reportPlots.MSPeakLists[grp] && reportPlots.MSPeakLists[grp].MSMS) || "";
        }
    }

    for (ann of [ "formulas", "compounds" ])
    {
        if (document.getElementById(ann + "Tab"))
        {
            Reactable.setFilter(ann + "Tab", "group", grp);
            
            // update susp only filter if needed
            const suspCheckEl = document.getElementById(ann + "-susp_only");
            if (suspCheckEl.style.display != "none" && suspCheckEl.checked)
            {
                // NOTE: we have to give the rowIndex below, as it's otherwise taken from the table's metadata,
                // which doesn't seem to be updated at this point...
                toggleAnnOnlySusp(ann, true, rowIndex);
            }

            if (ann === "compounds")
                document.getElementById('openMF').href = Reactable.getState('compoundsTab').meta.mfWebLinks[grp];
        }
    }
    
    const ccd = document.getElementById('comps_cluster-dendro');
    if (ccd)
    {
        ccd.src = reportPlots.compsCluster[grp].dendro;
        Array.from(document.getElementsByClassName('mcs')).forEach(el => el.style.display = (el.classList.contains('mcs-' + grp)) ? '' : 'none');
    }

    if (tabEl === "detailsTabSuspects")
    {
        const structEl = document.getElementById('struct_view-suspect');
        structEl.src = reportPlots.structs[rowValues.susp_InChIKey] || "";
        Reactable.setFilter('suspInfoTab', 'name', rowValues.susp_name);
        if (document.getElementById('suspAnnTab'))
            Reactable.setFilter('suspAnnTab', 'suspID', rowValues.susp_name + '-' + rowValues.group);
    }
    else if (tabEl === "detailsTabISTDs")
    {
        const structEl = document.getElementById('struct_view-istd');
        structEl.src = reportPlots.structs[rowValues.InChIKey] || "";
        Reactable.setFilter('ISTDInfoTab', 'name', rowValues.susp_name);
    }
    else if (tabEl === "detailsTabComponents")
    {
        let chromEl = document.getElementById('chrom_view-component');
        let specEl = document.getElementById('spectrum_view-component');
        let profileRelEl = document.getElementById('profileRel_view-component');
        let profileAbsEl = document.getElementById('profileAbs_view-component');
        const pl = reportPlots.components.components[rowValues.component];
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
        if (Object.keys(reportPlots.TPs).length > 0)
        {
            const specSimEl = document.getElementById('similarity_spec');
            specSimEl.src = reportPlots.TPs[rowValues.component][rowValues.cmpIndex - 1];
            specSimEl.style.display = ''; // may have been hidden if a previous img didn't exist
        }
        
        Reactable.setFilter('TPCandidatesTab', 'component', rowValues.component);
        Reactable.setFilter('TPCandidatesTab', 'group', rowValues.group);
        // activate first row
        // UNDONE: does this work properly with paging?
        const TPCandInstData = Array.from(Reactable.getInstance("TPCandidatesTab").data);
        const firstRowInd = TPCandInstData.findIndex(el => el.component === rowValues.component && el.group === rowValues.group);
        updateTPCandTabRowSel(TPCandInstData[firstRowInd], firstRowInd);
        
        if (document.getElementById('suspAnnTab'))
            Reactable.setFilter('suspAnnTab', 'suspID', rowValues.susp_name + '-' + rowValues.group);
        if (document.getElementById('similarityTab'))
            Reactable.setFilter('similarityTab', 'cmpID', rowValues.component + '-' + rowValues.group);
    }
    else if (tabEl === "detailsTabTPsParents")
    {
        document.getElementById("TPCompon-select").value = rowValues.component;
        updateTPCompon(rowValues.component, false);
    }
}

function updateTPCandTabRowSel(rowValues, rowIndex)
{
    Reactable.setMeta("TPCandidatesTab", { selectedRow: rowIndex });
}

function showTabCols(id, columnGroup, show)
{
    const cols = Reactable.getState(id).meta.colToggles[columnGroup];
    if (Array.isArray(cols))
        cols.forEach(col => Reactable.toggleHideColumn(id, col, !show));
    else
        Reactable.toggleHideColumn(id, cols, !show);
    
    // UNDONE: don't do this here
    if (column === "chrom_large")
        Reactable.toggleHideColumn(id, "chrom_small", show);
}

function updateTPCompon(cmpName, activateFG = true)
{
    const chromEl = document.getElementById('chrom_view-parent');
    const intEl = document.getElementById('int_plot-parent');
    
    if (chromEl)
        chromEl.src = reportPlots.chromsLarge[TPComponParentInfo[cmpName].group];
    if (intEl)
        intEl.src = reportPlots.intPlots[TPComponParentInfo[cmpName].group];
    
    if (document.getElementById('parentInfoTab'))
    {
        document.getElementById('struct_view-parent').src = reportPlots.structs[TPComponParentInfo[cmpName].InChIKey] || "";
        Reactable.setFilter('parentInfoTab', 'name', TPComponParentInfo[cmpName].name);
        // UNDONE: move to separate tab
        /*document.getElementById('struct_view-tp').src = reportPlots.structs[rowValues.susp_InChIKey] || "";
        Reactable.setFilter('TPInfoTab', 'name', rowValues.susp_name);*/
    }
    
    showTPGraph(cmpName);
    
    if (activateFG)
    {
        // activate proper FG table row: for the parents view switch to the selected component, otherwise to the first row
        const data = Array.from(Reactable.getInstance(getSelFGTableElement()).data);
        const firstRowInd = data.findIndex(el => el.component === cmpName);
        updateFeatTabRowSel(data[firstRowInd], firstRowInd);
    }
    
    if (!document.getElementById("viewTPDetailsParents").checked)
        Reactable.setFilter("detailsTabTPs", "component", cmpName);
}

function advanceTPCompon(dir)
{
    // based on https://stackoverflow.com/a/11556996
    let el = document.getElementById("TPCompon-select");
    let newIndex = el.selectedIndex + dir;
    if (newIndex < 0)
        newIndex = el.options.length - 1;
    else if (newIndex >= el.options.length)
        newIndex = 0;
    el.options[newIndex].selected = true;
    updateTPCompon(el.options[newIndex].value);
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

function filtColSelModalInit(tab, col)
{
    // UNDONE: this always adds all (unique) values, including those that are e.g. filtered out due to group selection    
    const vals = new Set(Reactable.getState(tab).data.map(row => row[col]));
    
    filtSelModalInit(tab, col, vals);
}

function filtFeatAnnSelModalInit(tab)
{
    filtSelModalInit(tab, "annotations", [ "None", "MS/MS", "Formulas", "Compounds" ]);
}

function filtSelModalInit(tab, col, vals)
{
    const mname = "filter_" + col;
    const curF = Reactable.getState(tab).meta[mname];

    let lg = document.getElementById("filterSelectListGroup");
    while (lg.firstChild)
        lg.removeChild(lg.firstChild);

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

function toggleConcsFilters(e)
{
    // as above, for feature table
    Reactable.getInstance("concsTab").allColumns.forEach(function(col)
    {
        if (col.name !== "group") // UNDONE: do this more elegantly?
            col.filterable = e;
    })
    
    applyFilterToggle("concsTab", "group", true);
}

function toggleToxFilters(e)
{
    // as above, for feature table
    Reactable.getInstance("toxTab").allColumns.forEach(function(col)
    {
        if (col.name !== "group") // UNDONE: do this more elegantly?
            col.filterable = e;
    })
    
    applyFilterToggle("toxTab", "group", true);
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

function toggleAnnOnlySusp(wh, e, r = undefined)
{
    const tid = (wh === "formulas") ? "formulasTab" : "compoundsTab";
    if (!e)
        Reactable.setFilter(tid, "suspect", undefined);
    else
    {
        const fgTab = getSelFGTableElement();
        const curRow = r || Reactable.getState(fgTab).meta.selectedRow;
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

function setTabGroupBy(tab, gb, expBtn = undefined)
{
    Reactable.setGroupBy(tab, [ gb ]);
    
    if (expBtn)
    {
        // hide toggle expand button if no grouping
        document.getElementById(expBtn).style.display = (gb === "") ? "none" : "";
    }
}

function downloadCSV(tab, out)
{
    const cols = Reactable.getState(tab).meta.CSVCols;
    Reactable.downloadDataCSV(tab, out, { columnIds: cols });
}


$(document).ready(function() {
    // Image zooming, based on https://stackoverflow.com/a/57694495
    $('body').prepend("<div class=\"zoomDiv\"><img src=\"\" class=\"zoomImg\"></div>");
    $('body').on('click', 'img:not(.zoomImg, .noZoomImg)', function() {
        const src = $(this).attr('data-srcZoom') || $(this).attr('src');
        $('.zoomImg').attr('src', src);
        $('.zoomDiv').css({opacity: '1', width: '70%'});
    });
    $('img.zoomImg').click(function() {
        $('.zoomDiv').css({opacity: '0', width: '0%'});
    });
    
    updateView("Plain");
    
    let el = document.getElementById("TPCompon-select")
    if (el)
        updateTPCompon(el.options[0].value, false);
    
    document.getElementById("filterRangeModal").addEventListener('shown.bs.modal', () => {
        document.getElementById("filtNumMin").focus();
    });
});
