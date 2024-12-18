// SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
//
// SPDX-License-Identifier: GPL-3.0-only

function selectPrevHot(which)
{
    // get rhot instance: https://github.com/jrowen/rhandsontable/issues/97
    var ht;
    if (which == "group")
        ht = HTMLWidgets.getInstance(groupHot).hot;
    else
        ht = HTMLWidgets.getInstance(metFragHot).hot;
    
    var nrow = ht.countRows();
    var r = Math.max(0, nrow - 1);
    if (ht.getSelected() && ht.getSelected()[0] > 0)
        r = ht.getSelected()[0] - 1;
    ht.selectCell(r, 0, undefined, undefined, undefined, true);
}

function selectNextHot(which)
{
    var ht;
    if (which == "group")
        ht = HTMLWidgets.getInstance(groupHot).hot;
    else
        ht = HTMLWidgets.getInstance(metFragHot).hot;

    var nrow = ht.countRows();
    var r = 0;
    if (ht.getSelected() && ht.getSelected()[0] < (nrow-1))
        r = ht.getSelected()[0] + 1;
    ht.selectCell(r, 0, undefined, undefined, undefined, true);
}
