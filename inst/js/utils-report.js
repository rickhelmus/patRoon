function getSpecPageElement()
{
    return document.getElementById("feature-groups").parentElement.parentElement;
}

function showSpec(group, type)
{
    clearComponentSpecs();
    
    var comps = document.getElementsByClassName('specClass');
    var specId = type + "_" + "spectra-" + group;
    var EICId = "EIC-" + group;
    for (var i=0; i<comps.length; i++)
        comps[i].style.display = (comps[i].classList.contains(specId) || comps[i].classList.contains(EICId)) ? 'flex' : 'none';

    document.getElementById("noSpecSelected").style.display = 'none';
}

function showComponentSpec(component, parentID)
{
    clearComponentSpecs();
    disableAllSpecs();
    
    var comp = document.getElementById(component).parentElement;
    var clonedComp = comp.cloneNode(true);
    clonedComp.id = clonedComp.id + "-clone";
    clonedComp.className += " clonedComp";
    
    getSpecPageElement().appendChild(clonedComp);
}

function clearComponentSpecs()
{
    $('.clonedComp').remove();
}

function disableAllSpecs()
{
    var comps = document.getElementsByClassName('specClass');
    for (var i=0; i<comps.length; i++)
        comps[i].style.display = 'none';
}

$(document).ready(function()
{
    //     setTimeout(disableAllSpecs, 500); // HACK: wait a bit so that HTML instances are available
});
