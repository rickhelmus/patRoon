// UNDONE: remove
function getGroupHot()
{
    // get rhot instance: https://github.com/jrowen/rhandsontable/issues/97
    return HTMLWidgets.getInstance(fGroupsHot).hot;
}

// UNDONE: remove
Shiny.addCustomMessageHandler("selectFGroupRow", function(row)
{
    getGroupHot().selectCell(row - 1, 0, undefined, undefined, undefined, true);
});

// UNDONE: remove
Shiny.addCustomMessageHandler("toggleFGroupRow", function(row)
{
    var ht = getGroupHot();
    ht.setDataAtCell(row - 1, 1, !ht.getDataAtCell(row - 1, 1));
});

function getPrimaryHot()
{
    // get rhot instance: https://github.com/jrowen/rhandsontable/issues/97
    return HTMLWidgets.getInstance(primaryHot).hot;
}

Shiny.addCustomMessageHandler("selectPrimaryRow", function(row)
{
    getPrimaryHot().selectCell(row - 1, 0, undefined, undefined, undefined, true);
});

Shiny.addCustomMessageHandler("togglePrimaryRow", function(row)
{
    var ht = getPrimaryHot();
    ht.setDataAtCell(row - 1, 1, !ht.getDataAtCell(row - 1, 1));
});

// show popup if something was changed and user is about to leave...
// based on https://stackoverflow.com/q/55048802
var sessionChanged = false;
Shiny.addCustomMessageHandler("setSessionChanged", function(changed) { sessionChanged = changed; });
window.addEventListener('beforeunload', function(e)
{
    // from https://stackoverflow.com/a/61404006
    if (sessionChanged)
    {
        e.preventDefault();
        e.returnValue = '';
        return;
    }

    delete e['returnValue'];
});
