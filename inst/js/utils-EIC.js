function getGroupHot()
{
    // get rhot instance: https://github.com/jrowen/rhandsontable/issues/97
    return HTMLWidgets.getInstance(fGroupHot).hot;
}

Shiny.addCustomMessageHandler("selectFGroupRow", function(row)
{
    getGroupHot().selectCell(row - 1, 0, undefined, undefined, undefined, true);
});

Shiny.addCustomMessageHandler("toggleFGroupRow", function(row)
{
    var ht = getGroupHot();
    ht.setDataAtCell(row - 1, 1, !ht.getDataAtCell(row - 1, 1));
});
