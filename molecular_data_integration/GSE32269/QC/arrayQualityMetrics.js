// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, true, false, true, false, false, true, false, false, false ];
var arrayMetadata    = [ [ "1", "GSE32269_1", "GSM799468_CL2003022549AA_215.CEL", "Tumour", "02/25/03 12:26:29" ], [ "2", "GSE32269_2", "GSM799469_CL2003022550AA_223.CEL", "Tumour", "02/25/03 12:39:19" ], [ "3", "GSE32269_3", "GSM799470_CL2003022558AA_358.CEL", "Tumour", "02/25/03 12:48:42" ], [ "4", "GSE32269_4", "GSM799471_CL2003022559AA_361.CEL", "Tumour", "02/25/03 13:00:17" ], [ "5", "GSE32269_5", "GSM799472_CL2003022560AA_368.CEL", "Tumour", "02/25/03 13:15:52" ], [ "6", "GSE32269_6", "GSM799473_CL2003022561AA_371.CEL", "Tumour", "02/25/03 13:39:11" ], [ "7", "GSE32269_7", "GSM799474_CL2003022551AA_224A.CEL", "Tumour", "02/25/03 12:49:59" ], [ "8", "GSE32269_8", "GSM799475_CL2003022552AA_224B.CEL", "Tumour", "02/25/03 11:40:06" ], [ "9", "GSE32269_9", "GSM799476_CL2003022553AA_228A.CEL", "Tumour", "02/25/03 13:04:01" ], [ "10", "GSE32269_10", "GSM799490_CL2004012114AA.CEL", "Metastasis", "01/21/04 15:51:27" ], [ "11", "GSE32269_11", "GSM799491_CL2004012124AA.CEL", "Metastasis", "01/21/04 14:41:00" ], [ "12", "GSE32269_12", "GSM799492_CL2004012125AA.CEL", "Metastasis", "01/21/04 15:40:10" ], [ "13", "GSE32269_13", "GSM799493_CL2004012127AA.CEL", "Metastasis", "01/21/04 15:29:02" ], [ "14", "GSE32269_14", "GSM799494_CL2004012128AA.CEL", "Metastasis", "01/21/04 15:16:33" ], [ "15", "GSE32269_15", "GSM799495_CL2004012131AA.CEL", "Metastasis", "01/21/04 15:05:52" ], [ "16", "GSE32269_16", "GSM799496_CL2004012136AA.CEL", "Metastasis", "01/21/04 14:40:09" ], [ "17", "GSE32269_17", "GSM799497_CL2004012139AA.CEL", "Metastasis", "01/21/04 14:29:20" ], [ "18", "GSE32269_18", "GSM799498_CL2004012145AA.CEL", "Metastasis", "01/21/04 16:39:06" ] ];
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
