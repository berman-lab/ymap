<?php
	session_start();
	require_once 'constants.php';
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html lang="en">
<head>
	<style type="text/css">
		html * {
			font-family: arial !important;
		}
		.upload {
			width:          100%;
			border:         0;
			height:         0px;
			vertical-align: middle;
			align:          left;
		}
		.tab {
			margin-left:    1cm;
		}
		.select {
			cursor: pointer;
		}
	</style>
	<meta http-equiv="content-type" content="text/html; charset=utf-8">
	<title>Y-MAP</title>
<!-- Used by secondary pages to update page on completion of processing. --!>
	<script type="text/javascript">
// Project/Dataset interface functions.
	function showColors(colorName,targetToChange,contentString) {
		if (document.getElementById(targetToChange)) {
			var select       = document.getElementById(targetToChange);
			select.innerHTML = '';
			fillerString     = '<div style=\"width:35px; height:20px; overflow:hidden; display:inline-block; background-color:';
			if        (colorName == "deep pink")       { fillerString = fillerString + 'rgb(255,0  ,127)';
			} else if (colorName == "magenta")         { fillerString = fillerString + 'rgb(255,0  ,255)';
			} else if (colorName == "electric indigo") { fillerString = fillerString + 'rgb(127,0  ,255)';
			} else if (colorName == "blue")            { fillerString = fillerString + 'rgb(0  ,0  ,255)';
			} else if (colorName == "dodger blue")     { fillerString = fillerString + 'rgb(0  ,127,255)';
			} else if (colorName == "cyan")            { fillerString = fillerString + 'rgb(0  ,255,255)';
			} else if (colorName == "spring green")    { fillerString = fillerString + 'rgb(0  ,255,127)';
			} else if (colorName == "green")           { fillerString = fillerString + 'rgb(0  ,255,0  )';
			} else if (colorName == "chartreuse")      { fillerString = fillerString + 'rgb(127,255,0  )';
			} else if (colorName == "yellow")          { fillerString = fillerString + 'rgb(255,255,0  )';
			} else if (colorName == "dark orange")     { fillerString = fillerString + 'rgb(255,127,0  )';
			} else if (colorName == "red")             { fillerString = fillerString + 'rgb(255,0  ,0  )';
			} else { fillerString = fillerString + 'rgb(127,127,127)';
			}
			fillerString     = fillerString + '\"><center><b>' + contentString + '</b></center></div>';
			select.innerHTML = fillerString;
		}
	}
	function update_project_remove_iframe(project_key, file_list) {
		project_key                          = project_key.replace('p_','');
		var show_button_element              = document.getElementById('panel_visualizeDataset_iframe').contentDocument.getElementById('show_'+project_key);
		show_button_element.setAttribute('data-file-list', file_list);
		show_button_element.style.visibility = 'visible';
		var project_iframe                   = document.getElementById('panel_manageDataset_iframe').contentDocument.getElementById(project_key);
		project_iframe.innerHTML             = '';
		document.getElementById('panel_manageDataset_iframe').contentDocument.getElementById("p_".project_key).style.height = "0px";
	}
	function update_project_label_color(project_key,label_color) {
		project_key                 = project_key.replace('p_','');
		var project_label1          = document.getElementById('panel_manageDataset_iframe').contentDocument.getElementById('p_label_'+project_key);
		var project_label2          = document.getElementById('panel_visualizeDataset_iframe').contentDocument.getElementById('p_label_'+project_key);
		project_label1.style.color  = label_color;
		project_label2.style.color  = label_color;
	}
	function update_project_file_size(project_key,sizeString_1,sizeString_2) {
		project_key                  = project_key.replace('p_','');
		var project_size1_span       = document.getElementById('panel_manageDataset_iframe').contentDocument.getElementById('p_size1_'+project_key);
		var project_size2_span       = document.getElementById('panel_manageDataset_iframe').contentDocument.getElementById('p_size2_'+project_key);
		if (sizeString_1 != "" && project_size1_span != null) {
			project_size1_span.innerHTML = " <font color='black' size='1'>("+sizeString_1+" bytes)</font>";
		}
		if (sizeString_2 != "" && project_size2_span != null) {
			project_size2_span.innerHTML = " <font color='black' size='1'>("+sizeString_2+" bytes)</font>";
		}
	}
	function resize_project(project_key, pixels) {
		document.getElementById('panel_manageDataset_iframe').contentDocument.getElementById(project_key).style.height = pixels+"px";
	}

// Hapmap interface functions.
	function update_hapmap_label_color(hapmap_key,label_color) {
		hapmap_key               = hapmap_key.replace('h_','');
		var hapmap_label         = document.getElementById('panel_hapmap_iframe').contentDocument.getElementById('h_label_'+hapmap_key);
		hapmap_label.style.color = label_color;
	}
	function resize_hapmap(hapmap_key, pixels) {
		document.getElementById('panel_hapmap_iframe').contentDocument.getElementById(hapmap_key).style.height = pixels+"px";
	}

// Genome interface functions.
	function update_genome_remove_iframe(genome_key) {
		console.log("$$ genome_key = '"+genome_key+"'");
		genome_key                           = genome_key.replace('g_','');
		var show_button_element              = document.getElementById('panel_genome_iframe').contentDocument.getElementById('show_'+genome_key);
		show_button_element.style.visibility = 'visible';
		var genome_iframe                    = document.getElementById('panel_genome_iframe').contentDocument.getElementById('frameContainer.g1_'+genome_key);
		genome_iframe.innerHTML              = '';
		document.getElementById("g_".genome_key).style.height = "0px";
	}
	function update_genome_label_color(genome_key,label_color) {
		genome_key               = genome_key.replace('g_','');
		var genome_label         = document.getElementById('panel_genome_iframe').contentDocument.getElementById('g_label_'+genome_key);
		genome_label.style.color = label_color;
	}
	function update_genome_file_size(genome_key,sizeString_1) {
		genome_key                  = genome_key.replace('g_','');
		var genome_size1_span       = document.getElementById('panel_genome_iframe').contentDocument.getElementById('g_size1_'+genome_key);
		if (sizeString_1 != "") {
			genome_size1_span.innerHTML = " <font color='black' size='1'>("+sizeString_1+" bytes)</font>";
		}
	}
	function resize_genome(genome_key, pixels) {
		document.getElementById('panel_genome_iframe').contentDocument.getElementById(genome_key).style.height = pixels+"px";
	}

<?php
	if (isset($_SESSION['logged_on'])) {
		$user = $_SESSION['user'];
	} else {
		$user = "default";
	}
?>
	user = "<?php echo $user ?>";
	</script>
</head>
<body>
<table width="100%"><tr>
	<td width="25%" align="center" style="max-height:100%">
		<table width="100%" height="300px"><tr valign="top"><td>
		<img src="images/Logo_title.2.png" alt="YMAP; Yeast Mapping Analysis Pipeline"><br><br>
		</td></tr><tr valign="bottom"><td align="middle">
<font size='2'>
	For support, please contact us at <a href="mailto:ymapsupport@tauex.tau.ac.il">ymapsupport@tauex.tau.ac.il</a><br/><br/>
<!-- Combined figures option is acting strangely. --!>
	<button onclick="Generate_combined_figure(); document.getElementById('combined_fig_options').style.display = 'inline';">Combine figures viewed below.</button><br>
<?php
	$cfig_CNV_SNP = "users/".$user."/combined_figure.1.png";
	$cfig_CNV     = "users/".$user."/combined_figure.2.png";
	$cfig_SNP     = "users/".$user."/combined_figure.3.png";

	$super_user_flag_file = "users/".$user."/super.txt";
	if (file_exists($super_user_flag_file)) {  // Super-user privilidges.
		$admin = "true";
	} else {
		$admin = "false";
	}
?>
	<div id='combined_fig_options' style='display:none;'>
		CNV-SNP/LOH <img src='images/icon_png_15b.png' alt-text='[PNG] button' align='center' onclick='loadExternal("<?php echo $cfig_CNV_SNP; ?>")'>
		CNV <img src='images/icon_png_15b.png' alt-text='[PNG] button' align='center' onclick='loadExternal("<?php echo $cfig_CNV; ?>")'>
		SNP/LOH <img src='images/icon_png_15b.png' alt-text='[PNG] button' align='center' onclick='loadExternal("<?php echo $cfig_SNP; ?>")'>
	</div>
</font>
		</td></tr></table>

<!--                                                 --!>
<!-- BEGIN : Construct tabbed menu and content pane. --!>
<!--                                                 --!>
</td><td width="75%" align="right" valign="top">
<table width="100%" height="<?php echo $ui_tabArea_height; ?>px" cellspacing="0">
<tr>
<td class="select" valign="middle" style="height:<?php echo $ui_tab_height; ?>; width:<?php echo $ui_tab_width; ?>;" align="center" id="tab_user"
	onclick="tabWindow('user');"            >User</td>
<td class="select" valign="middle" style="height:<?php echo $ui_tab_height; ?>; width:<?php echo $ui_tab_width; ?>;" align="center" id="tab_manageDataset"
	onclick="tabWindow('manageDataset');"   >Manage Datasets</td>
<td class="select" valign="middle" style="height:<?php echo $ui_tab_height; ?>; width:<?php echo $ui_tab_width; ?>;" align="center" id="tab_visualizeDataset"
	onclick="tabWindow('visualizeDataset');">Visualize Datasets</td>
<td class="select" valign="middle" style="height:<?php echo $ui_tab_height; ?>; width:<?php echo $ui_tab_width; ?>;" align="center" id="tab_genome"
	onclick="tabWindow('genome');"          >Reference Genome</td>
<td class="select" valign="middle" style="height:<?php echo $ui_tab_height; ?>; width:<?php echo $ui_tab_width; ?>;" align="center" id="tab_hapmap"
	onclick="tabWindow('hapmap');"          >Hapmap</td>
<td class="select" valign="middle" style="height:<?php echo $ui_tab_height; ?>; width:<?php echo $ui_tab_width; ?>;" align="center" id="tab_bugs"
	onclick="tabWindow('bugs');"            >Bug Reporting</td>
<td class="select" valign="middle" style="height:<?php echo $ui_tab_height; ?>; width:<?php echo $ui_tab_width; ?>;" align="center" id="tab_help"
	onclick="tabWindow('help');"            >Help</td>
<td class="select" valign="middle" style="height:<?php echo $ui_tab_height; ?>; width:<?php echo $ui_tab_width; ?>;" align="center" id="tab_examples"
	onclick="tabWindow('examples');"        >Example Datasets</td>
<?php
if (strcmp($admin,"true") == 0) {
	echo '<td class="select" valign="middle" style="height:'.$ui_tab_height.'; width:'.$ui_tab_width.';" align="center" id="tab_admin" ';
	echo 'onclick="tabWindow(\'admin\');"           >Admin</td>';
} else {
//	echo '<td class="select" valign="middle" style="height:'.$ui_tab_height.'; width:1px;" align="center" id="tab_admin">&nbsp;</td>';
	echo '<div id="tab_admin" style="visibility:hidden;"></div>';
}
?>

<td class="select" valign="middle" id="tab_blank">&nbsp;</td>
</tr>
<tr>
	<?php
		// $number_cols in the next line should be the number of tabs+1.
		$number_cols = 9;
		// If admin user, then one further column will be displayed.
		if (strcmp($admin,"true") == 0) {
			$number_cols += 1;
		}
	?>
	<td colspan="<?php echo $number_cols; ?>" valign="top" id="tab_content" name="tab_content" style="margin:0; padding:0;">
	<div id="panel_user"             name="panel_user"             style="margin:0; padding:0; border:none; width:100%; height:100%;"></div>
	<div id="panel_manageDataset"    name="panel_manageDataset"    style="margin:0; padding:0; border:none; width:100%; height:100%;"></div>
	<div id="panel_visualizeDataset" name="panel_visualizeDataset" style="margin:0; padding:0; border:none; width:100%; height:100%;"></div>
	<div id="panel_genome"           name="panel_genome"           style="margin:0; padding:0; border:none; width:100%; height:100%;"></div>
	<div id="panel_hapmap"           name="panel_hapmap"           style="margin:0; padding:0; border:none; width:100%; height:100%;"></div>
	<div id="panel_bugs"             name="panel_bugs"             style="margin:0; padding:0; border:none; width:100%; height:100%;"></div>
	<div id="panel_help"             name="panel_help"             style="margin:0; padding:0; border:none; width:100%; height:100%;"></div>
	<div id="panel_examples"         name="panel_examples"         style="margin:0; padding:0; border:none; width:100%; height:100%;"></div>
	<div id="panel_admin"            name="panel_admin"            style="margin:0; padding:0; border:none; width:100%; height:100%;"></div>
	</td>
</tr></table>
<!--                                               --!>
<!-- END : Construct tabbed menu and content pane. --!>
<!--                                               --!>
</td></tr></table>

<script type="text/javascript">
p_user                       = document.getElementById('panel_user');
p_manageDataset              = document.getElementById('panel_manageDataset');
p_visualizeDataset           = document.getElementById('panel_visualizeDataset');
p_genome                     = document.getElementById('panel_genome');
p_hapmap                     = document.getElementById('panel_hapmap');
p_bugs                       = document.getElementById('panel_bugs');
p_help                       = document.getElementById('panel_help');
p_examples                   = document.getElementById('panel_examples');
p_admin                      = document.getElementById('panel_admin');
p_user.innerHTML             = '<iframe id="panel_user_iframe"             src="panel.user.php"             style="margin:0; padding:0; border:none; width:100%; height:<?php echo $ui_iframe_height; ?>">';
p_manageDataset.innerHTML    = '<iframe id="panel_manageDataset_iframe"    src="panel.manageDataset.php"    style="margin:0; padding:0; border:none; width:100%; height:<?php echo $ui_iframe_height; ?>">';
p_visualizeDataset.innerHTML = '<iframe id="panel_visualizeDataset_iframe" src="panel.visualizeDataset.php" style="margin:0; padding:0; border:none; width:100%; height:<?php echo $ui_iframe_height; ?>">';
p_genome.innerHTML           = '<iframe id="panel_genome_iframe"           src="panel.genome.php"           style="margin:0; padding:0; border:none; width:100%; height:<?php echo $ui_iframe_height; ?>">';
p_hapmap.innerHTML           = '<iframe id="panel_hapmap_iframe"           src="panel.hapmap.php"           style="margin:0; padding:0; border:none; width:100%; height:<?php echo $ui_iframe_height; ?>">';
p_bugs.innerHTML             = '<iframe id="panel_bugs_iframe"             src="panel.bugs.php"             style="margin:0; padding:0; border:none; width:100%; height:<?php echo $ui_iframe_height; ?>">';
p_help.innerHTML             = '<iframe id="panel_help_iframe"             src="panel.help.php"             style="margin:0; padding:0; border:none; width:100%; height:<?php echo $ui_iframe_height; ?>">';
p_examples.innerHTML         = '<iframe id="panel_examples_iframe"         src="panel.examples.php"         style="margin:0; padding:0; border:none; width:100%; height:<?php echo $ui_iframe_height; ?>">';
p_admin.innerHTML            = '<iframe id="panel_admin_iframe"            src="panel.admin.php"            style="margin:0; padding:0; border:none; width:100%; height:<?php echo $ui_iframe_height; ?>;">';

function deselect_tab(name) {
	//console.log( '\tdeselect_tab("'+name+'")' );
	var current_tab = document.getElementById("tab_"+name);
	current_tab.style.border="1px solid #000000";
	if (name != "admin") {
		current_tab.style.backgroundColor="#DDDDDD";
	} else {
		current_tab.style.backgroundColor="#DDBBBB";
	}
	if (name != "user") {
		current_tab.style.borderLeft="none";
	}
	document.getElementById("panel_"+name).style.display='none';
}
function deselect_all_tabs() {
	//console.log( 'deselect_all_tabs()' );
	deselect_tab("user");
	deselect_tab("manageDataset");
	deselect_tab("visualizeDataset");
	deselect_tab("genome");
	deselect_tab("hapmap");
	deselect_tab("bugs");
	deselect_tab("help");
	deselect_tab("examples");
	deselect_tab("admin");
}
function select_tab(name) {
	//console.log ( 'select_tab("'+name+'")' );
	var current_tab = document.getElementById("tab_"+name);
	current_tab.style.border="1px solid #000000";
	current_tab.style.borderBottom="none";
	if (name != "admin") {
		current_tab.style.backgroundColor="#FFFFFF";
	} else {
		current_tab.style.backgroundColor="#FFDDDD";
	}
	if (name != "user") {
		current_tab.style.borderLeft="none";
	}
	document.getElementById("panel_"+name).style.display='inline';
}
function tabWindow(name) {
	deselect_all_tabs();
	select_tab(name);
	blank_and_content_tab();
	var tabInUse = name;
	localStorage.setItem("tabInUse", tabInUse);
}
function blank_and_content_tab() {
	var current_tab = document.getElementById("tab_blank");
	    current_tab.style.borderBottom="1px solid #000000";
	var current_tab = document.getElementById("tab_content");
	    current_tab.style.border="1px solid #000000";
	    current_tab.style.borderTop="none";
}


<!-- --!>
<!-- Project display section --!>
<!-- --!>
	function loadImage(key,imageUrl,imageScale) {
		document.getElementById('fig_'+key).innerHTML = "<img src='"+imageUrl+"' width='"+imageScale+"%'>";
	}
	function loadExternal(imageUrl) {
		window.open(imageUrl);
	}
	function showImg(imgSrc, H, W, Caption) {
		var newImg = window.open("","myImg",config="height="+H+",width="+W+"")
		newImg.document.write("<title>"+ Caption +"</title>")
		newImg.document.write("<img src='"+ imgSrc +"' height='"+ H +"' width='"+ W +"' onclick='window.close()' style='position:absolute;left:0;top:0'>")
		newImg.document.write("<script type='text/javascript'> document.oncontextmenu = new Function('return false') </script>")
		newImg.document.close();
	}
	function openProject(user,project,key,color1,color2,parent) {
		var visualize_iframe    = document.getElementById('panel_visualizeDataset_iframe');
		var show_button_element = visualize_iframe.contentDocument.getElementById("show_"+key);
		closeProject_viewOnly(key);
		console.log('#     parent.openProject : "'+user+':'+project+':'+key+':'+color1+':'+color2+':'+parent+'"');
		//console.log(show_button_element);

		if (show_button_element.checked == false) {
			closeProject(user,project,key,color1,color2,parent);
		} else {
			var file_list   = JSON.parse(show_button_element.getAttribute('data-file-list'));
			var file_prefix = "users/"+user+"/projects/"+project+"/";
			// Prefix all files with their folder path, so that we won't have to
			// manually add it in in every URL during HTML construction:
			for (fileIx = 0; fileIx < file_list.length; fileIx++) {
				file_list[fileIx] = file_prefix + file_list[fileIx];
			}

			var fig_linear_CNV_SNP               = file_prefix + "fig.CNV-SNP-map.2.";
			var fig_linear_CNV_SNP_RedGreen      = file_prefix + "fig.CNV-SNP-map.RedGreen.2.";
			var fig_standard_CNV_SNP             = file_prefix + "fig.CNV-SNP-map.1.";
			var fig_linear_CNV                   = file_prefix + "fig.CNV-map.2.";
			var fig_linear_CNV_highTop           = file_prefix + "fig.CNV-map.highTop.2.";
			var fig_standard_CNV                 = file_prefix + "fig.CNV-map.1.";
			if (file_list.indexOf(file_prefix + "fig.allelic_ratio-map.c2.png") != -1) {
				// ddRADseq.
				var fig_linear_SNP           = file_prefix + "fig.allelic_ratio-map.c2.";
				var fig_standard_SNP         = file_prefix + "fig.allelic_ratio-map.c1.";
			} else {
				// other.
				var fig_linear_SNP           = file_prefix + "fig.SNP-map.2.";
				var fig_standard_SNP         = file_prefix + "fig.SNP-map.1.";
			}
			var fig_linear_manual                = file_prefix + "fig.CNV-manualLOH-map.2.";
			var fig_standard_manual              = file_prefix + "fig.CNV-manualLOH-map.1.";
			var CGD_CNV_track                    = file_prefix + "cnv."+project+".gff3";
			var CGD_SNP_track                    = file_prefix + "allele_ratios."+project+".bed";

			var CNV_bias_SnpCghArray_GCcontent   = file_prefix + "fig_GCratio_vs_CGH.png";
			var CNV_bias_SnpCghArray_end         = file_prefix + "fig_EndDistance_vs_CGH.png";

			var CNV_bias_WGseq_end               = file_prefix + "fig.bias_chr_end.png";
			var CNV_bias_WGseq_GCcontent         = file_prefix + "fig.bias_GC_content.png";

			var CNV_bias_ddRADseq_1              = file_prefix + "fig.examine_bias.1.png";
			var CNV_bias_ddRADseq_2              = file_prefix + "fig.examine_bias.2.png";
			var CNV_bias_ddRADseq_3              = file_prefix + "fig.examine_bias.3.png";

			var fig_linear_SNPratio_histogram    = file_prefix + "fig.allelic_fraction_histogram.png";
			var fig_linear_SNPratio_fireplot     = file_prefix + "fig.allelic_ratio-map.b2.png";
			var visible_list                     = document.getElementById("visible_list");
			var string1 = "<div id='figure_"+key+"'><table border='0' align='center' width='100%'><tr><td width='35%' align='left'>";
			string1     = string1 + "<table><tr><td valign='top'>";
			string1     = string1 + project+" ";
			string1     = string1 + "</td><td valign='bottom'>";
			string1     = string1 + "<div id='userProjectA_"+key+"'   style='display:inline'></div>";
			string1     = string1 + "<div id='userProjectHET_"+key+"' style='display:inline'></div>";
			string1     = string1 + "<div id='userProjectB_"+key+"'   style='display:inline'></div>";
			string1     = string1 + "<div id='userProjectHOM_"+key+"' style='display:inline'></div>";
			string1     = string1 + "</td></tr></table>";
			string1     = string1 + "</td><td width='60%' align='left'><font size='-1'>";

			if (file_list.indexOf(fig_linear_CNV_SNP+"png") != -1) {
				mainFigure1 = fig_linear_CNV_SNP;
			} else if (file_list.indexOf(fig_linear_CNV+"png") != -1) {
				mainFigure1 = fig_linear_CNV;
			} else if (file_list.indexOf(fig_linear_SNP+"png") != -1) {
				mainFigure1 = fig_linear_SNP;
			}

			if (file_list.indexOf(fig_linear_CNV_SNP+"png") != -1) {
				string1 = string1 + "<b>CNV and SNP/LOH</b> (linear ";
				string1 = string1 + "<img src='images/icon_png_15b.png' alt-text='[PNG] button' align='center' onclick='loadImage(\""+key+"\",\""+fig_linear_CNV_SNP+"png\",\"100\")'> ";
				string1 = string1 + "<img src='images/icon_eps_15b.png' alt-text='[EPS] button' align='center' onclick='loadExternal(\""+fig_linear_CNV_SNP+"eps\")'>";
				if (file_list.indexOf(fig_linear_CNV_SNP_RedGreen+"png") != -1) {
					string1 = string1 + " or alternate colors ";
					string1 = string1 + "<img src='images/icon_png_15b.png' alt-text='[PNG] button' align='center' onclick='loadImage(\""+key+"\",\""+fig_linear_CNV_SNP_RedGreen+"png\",\"100\")'> ";
					string1 = string1 + "<img src='images/icon_eps_15b.png' alt-text='[EPS] button' align='center' onclick='loadExternal(\""+fig_linear_CNV_SNP_RedGreen+"eps\")'>";
				}
				string1 = string1 + " or standard ";
				string1 = string1 + "<img src='images/icon_png_15b.png' alt-text='[PNG] button' align='center' onclick='loadImage(\""+key+"\",\""+fig_standard_CNV_SNP+"png\",\"50\")'> ";
				string1 = string1 + "<img src='images/icon_eps_15b.png' alt-text='[EPS] button' align='center' onclick='loadExternal(\""+fig_standard_CNV_SNP+"eps\")'>";
				if (file_list.indexOf(fig_linear_manual+"png") != -1) {
					string1 = string1 + " or Linear-Manual ";
					string1 = string1 + "<img src='images/icon_png_15b.png' alt-text='[PNG] button' align='center' onclick='loadImage(\""+key+"\",\""+fig_linear_manual+"png\",\"100\")'> ";
					string1 = string1 + "<img src='images/icon_eps_15b.png' alt-text='[EPS] button' align='center' onclick='loadExternal(\""+fig_linear_manual+"eps\")'>";
				}
				string1 = string1 + ")";
				if ((file_list.indexOf(fig_linear_CNV+"png") != -1) || (file_list.indexOf(fig_linear_SNP+"png") != -1)) {
					string1 = string1 + " ";
				}
			}
			if (file_list.indexOf(fig_linear_CNV+"png") != -1) {
				string1 = string1 + "<br><b>CNV only</b> (lin. ";
				string1 = string1 + "<img src='images/icon_png_15b.png' alt-text='[PNG] button' align='center' onclick='loadImage(\""+key+"\",\""+fig_linear_CNV+"png\",\"100\")'> ";
				string1 = string1 + "<img src='images/icon_eps_15b.png' alt-text='[EPS] button' align='center' onclick='loadExternal(\""+fig_linear_CNV+"eps\")'>";
				string1 = string1 + " or sta. ";
				string1 = string1 + "<img src='images/icon_png_15b.png' alt-text='[PNG] button' align='center' onclick='loadImage(\""+key+"\",\""+fig_standard_CNV+"png\",\"50\")'> ";
				string1 = string1 + "<img src='images/icon_eps_15b.png' alt-text='[EPS] button' align='center' onclick='loadExternal(\""+fig_standard_CNV+"eps\")'>";
				string1 = string1 + ")";
				if (file_list.indexOf(fig_linear_SNP+"png") != -1) {
					string1 = string1 + " ";
				}
			}

			var CGD_tracks_present = false;

			// Show CNV bias figure for SnpCghArray, WGseq, and ddRADseq.
			if ((file_list.indexOf(CNV_bias_SnpCghArray_GCcontent) != -1) || (file_list.indexOf(CNV_bias_SnpCghArray_end) != -1)) {
				string1 = string1 + " : CNV biases ";
				if (file_list.indexOf(CNV_bias_SnpCghArray_GCcontent) != -1) {
					 string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+CNV_bias_SnpCghArray_GCcontent+"\",\"50\")'>%GC</button>";
				}
				if (file_list.indexOf(CNV_bias_SnpCghArray_end) != -1) {
					string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+CNV_bias_SnpCghArray_end+"\",\"50\")'>chr end</button>";
				}
			}
			if ((file_list.indexOf(CNV_bias_WGseq_end) != -1) || (file_list.indexOf(CNV_bias_WGseq_GCcontent) != -1)) {
				string1 = string1 + " : CNV biases ";
				if (file_list.indexOf(CNV_bias_WGseq_GCcontent) != -1) {
					string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+CNV_bias_WGseq_GCcontent+"\",\"50\")'>%GC</button>";
				}
				if (file_list.indexOf(CNV_bias_WGseq_end) != -1) {
					string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+CNV_bias_WGseq_end+"\",\"50\")'>chr end</button>";
				}
			}
			if (file_list.indexOf(fig_linear_CNV_highTop+"png") != -1) {
				string1 = string1 + " : CNV ";
				string1 = string1 +  "<button onclick='loadImage(\""+key+"\",\""+fig_linear_CNV_highTop+"png\",\"100\")'>high top</button>";
			}
			if ((file_list.indexOf(CNV_bias_ddRADseq_1) != -1) || (file_list.indexOf(CNV_bias_ddRADseq_2) != -1) || (file_list.indexOf(CNV_bias_ddRADseq_3) != -1)) {
				string1 = string1 + " : CNV biases ";
				if (file_list.indexOf(CNV_bias_ddRADseq_1) != -1) {
					string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+CNV_bias_ddRADseq_1+"\",\"100\")'>fragment length</button>";
				}
				if (file_list.indexOf(CNV_bias_ddRADseq_2) != -1) {
					string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+CNV_bias_ddRADseq_2+"\",\"100\")'>%GC</button>";
				}
				if (file_list.indexOf(CNV_bias_ddRADseq_3) != -1) {
					string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+CNV_bias_ddRADseq_3+"\",\"100\")'>chr end</button>";
				}
			}
			if (file_list.indexOf(CGD_CNV_track) != -1) {
				string1 += "<a href=\"" + CGD_CNV_track + "\">GBrowse CNV track</a>";
				CGD_tracks_present = true;
			}
			if (file_list.indexOf(fig_linear_SNP+"png") != -1) {
				string1 = string1 + "<br><b>SNP/LOH only</b> (lin. ";
				string1 = string1 + "<img src='images/icon_png_15b.png' alt-text='[PNG] button' align='center' onclick='loadImage(\""+key+"\",\""+fig_linear_SNP+"png\",\"100\")'> ";
				string1 = string1 + "<img src='images/icon_eps_15b.png' alt-text='[EPS] button' align='center' onclick='loadExternal(\""+fig_linear_SNP+"eps\")'>";
				string1 = string1 + " or sta. ";
				string1 = string1 + "<img src='images/icon_png_15b.png' alt-text='[PNG] button' align='center' onclick='loadImage(\""+key+"\",\""+fig_standard_SNP+"png\",\"50\")'> ";
				string1 = string1 + "<img src='images/icon_eps_15b.png' alt-text='[EPS] button' align='center' onclick='loadExternal(\""+fig_standard_SNP+"eps\")'>";
				string1 = string1 + ")";
			}

			// Show allelic ratio plot version for ddRADseq and WGseq.
			if ((file_list.indexOf(fig_linear_SNPratio_histogram) != -1) || (file_list.indexOf(fig_linear_SNPratio_fireplot) != -1)) {
				string1 = string1 + " : SNP ratios ";
				if (file_list.indexOf(fig_linear_SNPratio_histogram) != -1) {
					string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+fig_linear_SNPratio_histogram+"\",\"100\")'>histogram</button>";
				}
				if (file_list.indexOf(fig_linear_SNPratio_fireplot) != -1) {
					string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+fig_linear_SNPratio_fireplot+"\",\"100\")'>fire plot</button>";
				}
			}

			if (file_list.indexOf(CGD_SNP_track) != -1) {
				string1 += "<a href=\"" + CGD_SNP_track + "\">GBrowse allele ratio track</a>";
				CGD_tracks_present = true;
			}

			if (CGD_tracks_present)
			{
				string1 += "<br><a href=\"https://github.com/berman-lab/ymap/wiki/How-to-import-Ymap-generated-tracks-into-CGD's-GBrowse\">(help on importing tracks into GBrowse)</a>";
			}

			string1 = string1 + "</font>";

			var string2 = "</td><td width='5%' align='right'><div onclick=\'closeProject(\""+user+"\",\""+project+"\",\""+key+"\",\""+color1+"\",\""+color2+"\",\""+parent+"\");' style='display:inline-block;'><b>[X]</b></div></td></tr>";
			var string3 = "<tr><td align='center' colspan='3'>";
			var string4 = "<div id='fig_"+key+"'><img src='"+mainFigure1+"png' width='100%'></div>";
			string4 = string4 + "</td></tr></table>"+"<hr></div>";
			visible_list.innerHTML += string1+string2+string3+string4;
			if (color1 != 'null') { // Hapmap analysis done.
				showColors(color1 ,'userProjectA_'+key,'a');
				showColors(color2 ,'userProjectB_'+key,'b');
				showColors('red' ,'userProjectHOM_'+key,'hom');
				showColors('grey','userProjectHET_'+key,'ab');
			} else { // No hapmap analysis done.
				if (parent == project) { // no LOH calculations done.
					showColors('grey','userProjectHET_'+key,'ab');
				} else { // LOH calculations have been done vs. parent dataset.
					showColors('red' ,'userProjectHOM_'+key,'hom');
					showColors('grey','userProjectHET_'+key,'ab');
				}
			}
			var projectsShown = localStorage.getItem("projectsShown");
			//console.log("'"+projectsShown+"'");
			if (projectsShown != null) {
				projectsShown = projectsShown.replace(user+":"+project+":"+key+":null:null:null");
			} else {
				projectsShown = "";
			}
			projectsShown = projectsShown+" "+user+":"+project+":"+key+":null:null:null";
			projectsShown = projectsShown.replace("  "," ");   // remove duplicate " " characters.
			while (projectsShown.charAt(0) == " ")
				projectsShown = projectsShown.slice( 1 );      // remove leading " " character.
			localStorage.setItem("projectsShown", projectsShown);
			console.log('#     Add to projectsShown : "'+user+':'+project+':'+key+':null:null:null"');
			console.log('#         projectsShown = "'+projectsShown+'"');
		}
	}
	function closeProject(user,project,key,color1,color2,parent) {
		var visualize_iframe    = document.getElementById('panel_visualizeDataset_iframe');
        var show_button_element = visualize_iframe.contentDocument.getElementById("show_"+key);
		if (show_button_element) {
			show_button_element.checked = false;
		}
		var figure_element = document.getElementById("figure_"+key);
		if (figure_element) {
			figure_element.remove();
		}
		var projectsShown = localStorage.getItem("projectsShown");
		projectsShown = projectsShown.replace(user+":"+project+":"+key+":null:null:null","");
		projectsShown = projectsShown.replace("  "," ");  // remove duplicate " " characters.
		while (projectsShown.charAt(0) == " ")
			projectsShown = projectsShown.slice( 1 );     // remove leading " " characater.
		localStorage.setItem("projectsShown", projectsShown);
		console.log('#     Remove from projectsShown : "'+user+':'+project+':'+key+':null:null:null"');
		console.log('#         projectsShown = "'+projectsShown+'"');
	}
	function closeProject_viewOnly(key) {
		var figure_element = document.getElementById("figure_"+key);
		if (figure_element) {
			figure_element.remove();
		}
	}
	function CloseAllProjects(user) {
		var visualize_iframe    = document.getElementById('panel_visualizeDataset_iframe');
		// Function not implementd.
	}
	</script>
<hr>
<div id="visible_list" name="visible_list"></div>


<!=====================================================================!>
<?php if (isset($_SESSION['logged_on'])) { ?>
<DIV id="Hidden_InstallNewDataset" style="display: none; position: absolute; border: solid black 1px; padding: 5px; text-align: justify;">
	<table width="100%"><tr>
	<td width="95%" align="left">Install New Dataset</td>
	<td width="5%" align="right"><div onmousedown="document.getElementById('Hidden_InstallNewDataset').style.display = 'none';" style="display:inline-block;"><b>[X]</b></div></td>
	</tr></table>
	<br>
	<iframe id="Hidden_InstallNewDataset_Frame" src="project.create_window.php"></iframe>
</DIV>
<DIV id="Hidden_InstallNewGenome" style="display: none; position: absolute; border: solid black 1px; padding: 5px; text-align: justify;">
	<table width="100%"><tr>
	<td width="95%" align="left">Install New Genome</td>
	<td width="5%" align="right"><div onmousedown="document.getElementById('Hidden_InstallNewGenome').style.display = 'none';" style="display:inline-block;"><b>[X]</b></div></td>
	</tr></table>
	<br>
	<iframe id="Hidden_InstallNewGenome_Frame" src="genome.create_window.php"></iframe>
</DIV>
<DIV id="Hidden_GenerateNewHapmap" style="display: none; position: absolute; border: solid black 1px; padding: 5px; text-align: justify;">
	<table width="100%"><tr>
	<td width="95%" align="left">Generate New Hapmap</td>
	<td width="5%" align="right"><div onmousedown="document.getElementById('Hidden_GenerateNewHapmap').style.display = 'none';" style="display:inline-block;"><b>[X]</b></div></td>
	</tr></table>
	<br>
	<iframe id="Hidden_GenerateNewHapmap_Frame" src="hapmap.create_window.php"></iframe>
</DIV>
<DIV id="Hidden_AddToHapmap" style="display: none; position: absolute; border: solid black 1px; padding: 5px; text-align: justify;">
	<table width="100%"><tr>
	<td width="95%" align="left">Generate New Hapmap</td>
	<td width="5%" align="right"><div onmousedown="document.getElementById('Hidden_GenerateNewHapmap').style.display = 'none';" style="display:inline-block;"><b>[X]</b></div></td>
	</tr></table>
	<br>
	<iframe id="Hidden_AddToHapmap_Frame"></iframe>
</DIV>
<DIV id="Hidden_BugTracker" style="display: none; position: absolute; border: solid black 1px; padding: 5px; text-align: justify;">
	<table width="100%"><tr>
	<td width="95%" align="left">Bug Reports</td>
	<td width="5%" align="right"><div onmousedown="document.getElementById('Hidden_BugTracker').style.display = 'none';" style="display:inline-block;"><b>[X]</b></div></td>
	</tr></table>
	<br>
	<iframe id="Hidden_BugTracker_Frame" src="bug_window.php"></iframe>
</DIV>
<?php } ?>

<script type="text/javascript">
function show_hidden(hidden_panel) {
	var el,x,y,w,h, Fel,Fx,Fy,Fw,Fh;
	el = document.getElementById(hidden_panel);
	x  = 50;
	y  = 50;
	w  = window.innerWidth-100;
	h  = window.innerHeight-100;
	el.style.left            = x + "px";
	el.style.top             = y + "px";
	el.style.display         = "block";
	el.style.width           = w + "px";
	el.style.height          = h + "px";
	el.style.backgroundColor = "rgb(200,200,200)";
	Fel = document.getElementById(hidden_panel+"_Frame");
	Fx  = 25;
	Fy  = 25;
	Fw  = w-50;
	Fh  = h-50;
	Fel.style.left            = Fx + "px";
	Fel.style.top             = Fy + "px";
	Fel.style.width           = Fw + "px";
	Fel.style.height          = Fh + "px";
	Fel.style.position        = "absolute";
	Fel.style.border          = "0 px";
	Fel.style.backgroundColor = "rgb(255,255,255)";
}
function reload_hidden(hidden_panel,source) {
	Fel     = document.getElementById(hidden_panel+"_Frame");
	Fel.src = source;
}

// =====================================================================

function update_interface_after_login() {
	document.getElementById('panel_manageDataset_iframe'   ).contentWindow.location.reload();
	document.getElementById('panel_visualizeDataset_iframe').contentWindow.location.reload();
	document.getElementById('panel_genome_iframe'          ).contentWindow.location.reload();
	document.getElementById('panel_hapmap_iframe'          ).contentWindow.location.reload();
	document.getElementById('panel_bugs_iframe'            ).contentWindow.location.reload();
	document.getElementById('panel_help_iframe'            ).contentWindow.location.reload();
	document.getElementById('panel_examples_iframe'        ).contentWindow.location.reload();
}
function update_interface_after_logout() {
	localStorage.removeItem('projectsShown');
	var ff1 = parent.document.getElementById('panel_manageDataset_iframe');
	var ff2 = parent.document.getElementById('panel_visualizeDataset_iframe');
	var ff3 = parent.document.getElementById('panel_genome_iframe');
	var ff4 = parent.document.getElementById('panel_hapmap_iframe');
	var ff5 = parent.document.getElementById('panel_bugs_iframe');
	var ff6 = parent.document.getElementById('panel_help_iframe');
	var ff7 = parent.document.getElementById('panel_examples_iframe');
	ff1.src = ff1.src;  // reload panel.
	ff2.src = ff2.src;
	ff3.src = ff3.src;
	ff4.src = ff4.src;
	ff5.src = ff5.src;
	ff6.src = ff6.src;
	ff7.src = ff7.src;
}
function update_projectsShown_after_new_project() {
	projectsShown = localStorage.getItem("projectsShown");
	console.log('## before projectsShown = "'+projectsShown+'"');
	var projectsShown_entries = projectsShown.split(' ');
	var new_projectsShown = "";
	localStorage.setItem("projectsShown","");
	for (var i=0;i<projectsShown_entries.length; i++) {
		var currentProject = projectsShown_entries[i];
		if (currentProject != '') {
			var entry_parts    = currentProject.split(':');
			projID = parseInt(entry_parts[2])+1;
			projID.toString();
			new_project = entry_parts[0]+":"+entry_parts[1]+":"+projID+":"+entry_parts[3]+":"+entry_parts[4]+":"+entry_parts[5];
			console.log('## old_project = "'+currentProject+'"');
			console.log('## new_project = "'+new_project+'"');
			projectsShown = projectsShown.replace(currentProject, new_project);
		}
	}
	console.log('## after projectsShown = "'+projectsShown+'"');
	localStorage.setItem("projectsShown",projectsShown);
}
function update_projectsShown_after_project_delete(deletedProjectKey) {
	projectsShown = localStorage.getItem("projectsShown");
	// Clean up 'deletedProjectID' to projectID.   'p1_[number]_delete' => [number]
	var deletedProjectID = deletedProjectKey;
	    deletedProjectID = deletedProjectID.replace("p1_","");
	    deletedProjectID = deletedProjectID.replace("_delete","");
	    deletedProjectID = Number(deletedProjectID);
	// adjust projectsShown string to reflect new positions of projects.
	var projectsShown_entries = projectsShown.split(' ');
	var new_projectsShown = "";
	localStorage.setItem("projectsShown","");
	for (var i=0;i<projectsShown_entries.length; i++) {
		var currentProject = projectsShown_entries[i];
		if (currentProject != '') {
			var entry_parts    = currentProject.split(':');
			projID = parseInt(entry_parts[2]);
			if (projID < deletedProjectID) {
				projID.toString();
				new_projectsShown = new_projectsShown + entry_parts[0]+":"+entry_parts[1]+":"+projID+":"+entry_parts[3]+":"+entry_parts[4]+":"+entry_parts[5]+" ";
			} else if (projID == deletedProjectID) {
				Deleted_entry_parts = entry_parts;
				// Remove deleted entry from active area.
				closeProject(Deleted_entry_parts[0],Deleted_entry_parts[1],Deleted_entry_parts[2],Deleted_entry_parts[3],Deleted_entry_parts[4],Deleted_entry_parts[5]);
			} else {  // if (projID > deletedProjectID) {
				projID = projID-1;
				projID.toString();
				new_projectsShown = new_projectsShown + entry_parts[0]+":"+entry_parts[1]+":"+projID+":"+entry_parts[3]+":"+entry_parts[4]+":"+entry_parts[5]+" ";
			}
		}
	}
	new_projectsShown = new_projectsShown.slice(0, -1);
	localStorage.setItem("projectsShown",new_projectsShown);
}


// upon page reload.
if(localStorage.getItem("tabInUse")){
	var tabInUse      = localStorage.getItem("tabInUse");
} else {
	var tabInUse      = "user";
}
if(localStorage.getItem("projectsShown")){
	var projectsShown = localStorage.getItem("projectsShown");
}

// Reload previously viewed tab after page reload.
<?php
	if (!isset($_SESSION['logged_on'])) {
		echo "tabInUse = 'user';\n";
		echo "key_offset = '';\n";
	}
?>
console.log('#: tabInUse = "'+tabInUse+'"');
tabWindow(tabInUse);

// Reload previously viewed project datasets after page reload.
console.log('#: projectsShown = "'+projectsShown+'"');
function restore_shown_figures() {
	// function is called at the end of page loading.
	if (projectsShown) {
		var projectsShown_entries = projectsShown.split(' ');
		localStorage.setItem("projectsShown","");

		var visualize_iframe = document.getElementById('panel_visualizeDataset_iframe');
		visualize_iframe.onload = function () {
			console.log("#: Opening projects shown before page reload.");
			for (var i=0;i<projectsShown_entries.length; i++) {
				var currentProject = projectsShown_entries[i];
				console.log('#:    currentProject = '+currentProject);
				if (currentProject != '') {
					var entry_parts    = currentProject.split(':');

					// select checkBoxes for previously viewed datasets.
					key = entry_parts[2];
					//console.log('#:    4 elementID = show_'+key);

					var show_button_element = visualize_iframe.contentDocument.getElementById("show_"+key);
					//console.log(show_button_element);
					show_button_element.checked = true;

					// Open projects previously shown.
					//console.log('#:5 entry to show = "'+currentProject+'"');
					var entry_parts    = currentProject.split(':');
					openProject(entry_parts[0], entry_parts[1], entry_parts[2], 'null', 'null', 'null');
				}
			}
		}
	}
}

// ====== Page reload when logged out =================================

function reload_page() {
	var go = false;
	// trigger reload if reload_page is called and user is found to be logged off.
	<?php
	if (!isset($_SESSION['logged_on'])) {
		echo "\t// User is not logged in.\n";
		//location.reload();
	} else {
		echo "// User is still logged in.\n";
	}
	?>
}

<?php
if (isset($_SESSION['logged_on'])) {
	echo "\n// Initiate recurrent call to reload_page function, which depends upon project status.\n";
	echo "var intervalID = window.setInterval(reload_page, 3000);\n";
}
?>
</script>
<!-- Manage interfaces for deleting objects.--!>
	<script type="text/javascript" src="js/project.delete_index.js"></script>
	<script type="text/javascript" src="js/genome.delete_index.js"></script>
	<script type="text/javascript" src="js/user.delete_index.js"></script>
	<script type="text/javascript" src="js/hapmap.delete_index.js"></script>
<!-- Needed for deletion functions.--!>
	<script type="text/javascript" src="js/jquery-1.7.2.js"></script>
	<script type="text/javascript" src="js/jquery.form.js"></script>


<!-- Setup for Generating combined figures. --!>
<style>
	.hide { position:absolute; top:-1px; left:-1px; width:1px; height:1px; }
</style>
<iframe name="hiddenFrame" class="hide"></iframe>
<script type="text/javascript">
function Generate_combined_figure() {
	var user           = '<?php echo $user; ?>';
	var projectsShown  = localStorage.getItem("projectsShown");
	var autoSubmitForm = document.createElement("form");
	    autoSubmitForm.setAttribute("method","post");
	    autoSubmitForm.setAttribute("action","images.combine.php");
	    autoSubmitForm.setAttribute("target","hiddenFrame");
	var input2         = document.createElement("input");
	    input2.setAttribute("type","hidden");
	    input2.setAttribute("name","projectsShown");
	    input2.setAttribute("value",projectsShown);
	    autoSubmitForm.appendChild(input2);
	document.body.appendChild(autoSubmitForm);
	autoSubmitForm.submit();
}
restore_shown_figures();

function hide_combined_fig_menu() {
    menu = document.getElementById("combined_fig_options");
    menu.style.display = "none";
}
</script>



</body>
</html>
