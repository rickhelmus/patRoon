Shiny.addCustomMessageHandler("selectFGroupRow", function(row)
{
    // get rhot instance: https://github.com/jrowen/rhandsontable/issues/97
    var ht = HTMLWidgets.getInstance(groupHot).hot;
    ht.selectCell(row - 1, 0, undefined, undefined, undefined, true);
});
