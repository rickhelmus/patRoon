function selectPrevFGroup()
{
    // get rhot instance: https://github.com/jrowen/rhandsontable/issues/97
    var ht = HTMLWidgets.getInstance(groupHot).hot;
    var nrow = ht.countRows();
    var r = Math.max(0, nrow - 1);
    if (ht.getSelected() && ht.getSelected()[0][0] > 0)
        r = ht.getSelected()[0][0] - 1;
    ht.selectCell(r, 0, undefined, undefined, undefined, true);
}

function selectNextFGroup()
{
    var ht = HTMLWidgets.getInstance(groupHot).hot;
    var nrow = ht.countRows();
    var r = 0;
    if (ht.getSelected() && ht.getSelected()[0][0] < (nrow-1))
        r = ht.getSelected()[0][0] + 1;
    ht.selectCell(r, 0, undefined, undefined, undefined, true);
}

function setGroupHotContextMenu()
{
    var ht = HTMLWidgets.getInstance(groupHot).hot;
    ht.updateSettings(
    {
        contextMenu:
        {
            items:
            {
                "enableAll":
                {
                    name: "Enable all",
                    callback: function(key, options) { Shiny.onInputChange("enableAllGroups", Math.random()); }
                },
                "disableAll":
                {
                    name: "Disable all",
                    callback: function(key, options) { Shiny.onInputChange("disableAllGroups", Math.random()); }
                },
                
                "hsep1": "---------",
                
                "addDAEIC":
                {
                    key: "addDAEIC",
                    name: "Add DA EIC",
                    submenu:
                    {
                        items:
                        [
                            {
                                name: "to enabled analyses",
                                callback: function(key, options) { Shiny.onInputChange("addDAEICs", ["enabled", 0, Math.random()]); },
                                key: "addDAEIC:1"
                            }, 
                            {
                                name: "to all analyses",
                                callback: function(key, options) { Shiny.onInputChange("addDAEICs", ["all", 0, Math.random()]); },
                                key: "addDAEIC:2"
                            }, 
                        ]
                    }
                },
                "addDAEICBG":
                {
                    key: "addDAEICBG",
                    name: "Add DA EIC with bg subtr",
                    submenu:
                    {
                        items:
                        [
                            {
                                name: "to enabled analyses",
                                callback: function(key, options) { Shiny.onInputChange("addDAEICs", ["enabled", 1, Math.random()]); },
                                key: "addDAEICBG:1"
                            }, 
                            {
                                name: "to all analyses",
                                callback: function(key, options) { Shiny.onInputChange("addDAEICs", ["all", 1, Math.random()]); },
                                key: "addDAEICBG:2"
                            }, 
                        ]
                    }
                }
            }
        }
    })
}

function setAnalysesHotContextMenu()
{
    var ht = HTMLWidgets.getInstance(analysesHot).hot;
    ht.updateSettings(
    {
        contextMenu:
        {
            items:
            {
                "enableAll":
                {
                    name: "Enable all",
                    callback: function(key, options) { Shiny.onInputChange("enableAllAnalyses", Math.random()); }
                },
                "disableAll":
                {
                    name: "Disable all",
                    callback: function(key, options) { Shiny.onInputChange("disableAllAnalyses", Math.random()); }
                },
                
                "hsep1": "---------",
                
                "addDAEIC":
                {
                    key: "addDAEIC",
                    name: "Add DA EIC",
                    submenu:
                    {
                        items:
                        [
                            {
                                name: "to selected analysis",
                                callback: function(key, options) { Shiny.onInputChange("addDAEICs", ["selected", 0, Math.random()]); },
                                key: "addDAEIC:1"
                            }, 
                            {
                                name: "to enabled analyses",
                                callback: function(key, options) { Shiny.onInputChange("addDAEICs", ["enabled", 0, Math.random()]); },
                                key: "addDAEIC:1"
                            }, 
                            {
                                name: "to all analyses",
                                callback: function(key, options) { Shiny.onInputChange("addDAEICs", ["all", 0, Math.random()]); },
                                key: "addDAEIC:2"
                            },
                        ]
                    }
                },
                "addDAEICBG":
                {
                    key: "addDAEICBG",
                    name: "Add DA EIC with bg subtr",
                    submenu:
                    {
                        items:
                        [
                            {
                                name: "to selected analysis",
                                callback: function(key, options) { Shiny.onInputChange("addDAEICs", ["selected", 1, Math.random()]); },
                                key: "addDAEIC:1"
                            },                         
                            {
                                name: "to enabled analyses",
                                callback: function(key, options) { Shiny.onInputChange("addDAEICs", ["enabled", 1, Math.random()]); },
                                key: "addDAEICBG:1"
                            }, 
                            {
                                name: "to all analyses",
                                callback: function(key, options) { Shiny.onInputChange("addDAEICs", ["all", 1, Math.random()]); },
                                key: "addDAEICBG:2"
                            }, 
                        ]
                    }
                }
            }
        }
    })
}

function init()
{
    setGroupHotContextMenu()
    setAnalysesHotContextMenu()
    
    // HACK: use timeouts to make sure context menu is changed after everything has settled
    $('#groupHot').on('shiny:value', function(event) { setTimeout(setGroupHotContextMenu, 500); });
    $('#analysesHot').on('shiny:value', function(event) { setTimeout(setAnalysesHotContextMenu, 500); });
}

$(document).ready(function()
{
    // UNDONE: disable for now
    //setTimeout(init, 500); // HACK: wait a bit so that HTML instances are available
});
