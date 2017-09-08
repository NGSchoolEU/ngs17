// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, true, false, false, false, false ];
var arrayMetadata    = [ [ "1", "GSE8218_1", "GSM203263.CEL", "Normal", "02/06/03 13:24:06" ], [ "2", "GSE8218_2", "GSM203264.CEL", "Normal", "02/06/03 13:34:25" ], [ "3", "GSE8218_3", "GSM203266.CEL", "Normal", "02/06/03 13:55:12" ], [ "4", "GSE8218_4", "GSM203268.CEL", "Normal", "02/06/03 14:16:00" ], [ "5", "GSE8218_5", "GSM203269.CEL", "Normal", "02/06/03 14:37:31" ], [ "6", "GSE8218_6", "GSM203270.CEL", "Normal", "02/06/03 14:49:06" ], [ "7", "GSE8218_7", "GSM203271.CEL", "Normal", "02/06/03 14:59:17" ], [ "8", "GSE8218_8", "GSM203272.CEL", "Normal", "02/06/03 15:09:26" ], [ "9", "GSE8218_9", "GSM203276.CEL", "Normal", "02/04/03 13:12:46" ], [ "10", "GSE8218_10", "GSM203277.CEL", "Normal", "02/04/03 13:23:05" ], [ "11", "GSE8218_11", "GSM203309.CEL", "Tumour", "09/10/02 13:31:30" ], [ "12", "GSE8218_12", "GSM203314.CEL", "Tumour", "04/23/03 13:49:31" ], [ "13", "GSE8218_13", "GSM203324.CEL", "Tumour", "04/29/03 12:44:22" ], [ "14", "GSE8218_14", "GSM203333.CEL", "Tumour", "05/20/03 16:35:04" ], [ "15", "GSE8218_15", "GSM203334.CEL", "Tumour", "05/20/03 16:24:14" ], [ "16", "GSE8218_16", "GSM203338.CEL", "Tumour", "05/20/03 15:43:23" ], [ "17", "GSE8218_17", "GSM203341.CEL", "Tumour", "05/20/03 16:55:50" ], [ "18", "GSE8218_18", "GSM203396.CEL", "Tumour", "12/11/03 13:48:38" ], [ "19", "GSE8218_19", "GSM203398.CEL", "Tumour", "12/11/03 14:08:54" ], [ "20", "GSE8218_20", "GSM203400.CEL", "Tumour", "12/11/03 14:29:45" ] ];
var svgObjectNames   = [ "pca", "dens", "dig" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
	var success = false;
	i = 0; 
	/* Some of this looping could already be cached in reportInit() */
	while( (!success) & (i < ssrules.length) ) {
	    selector = ssrules[i].selectorText;  // The selector 
            if (!selector) 
		continue; // Skip @import and other nonstyle rules
            if (selector == (".aqm" + reportObjId)) {
		success = true; 
		ssrules[i].style.cssText = cssText[0+status];
	    } else {
		i++;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
