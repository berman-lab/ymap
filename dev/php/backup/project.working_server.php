<?php
    session_start();
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
		"http://www.w3.org/TR/html4/loose.dtd">
<head>
<style type="text/css">
body {font-family: arial;}
	.tab {
		margin-left:    1cm;
	}
	.clock {
		float:          left;
		margin-right:   0.25cm;
	}
</style>
</head>
<?php
	require_once 'constants.php';

	$user     = filter_input(INPUT_POST, "user",     FILTER_SANITIZE_STRING);
	$project  = filter_input(INPUT_POST, "project",  FILTER_SANITIZE_STRING);
	$key      = filter_input(INPUT_POST, "key",      FILTER_SANITIZE_STRING);
	$status   = filter_input(INPUT_POST, "status",   FILTER_SANITIZE_STRING);

	// increment clock animation...
	$status   = ($status + 1) % 12;
	if ($status == 0) {          $clock = "<img src=\"../images/12.png\" alt-text=\"12\" class=\"clock\" >";
	} else if ($status == 1) {   $clock = "<img src=\"../images/01.png\" alt-text=\"01\" class=\"clock\" >";
	} else if ($status == 2) {   $clock = "<img src=\"../images/02.png\" alt-text=\"02\" class=\"clock\" >";
	} else if ($status == 3) {   $clock = "<img src=\"../images/03.png\" alt-text=\"03\" class=\"clock\" >";
	} else if ($status == 4) {   $clock = "<img src=\"../images/04.png\" alt-text=\"04\" class=\"clock\" >";
	} else if ($status == 5) {   $clock = "<img src=\"../images/05.png\" alt-text=\"05\" class=\"clock\" >";
	} else if ($status == 6) {   $clock = "<img src=\"../images/06.png\" alt-text=\"06\" class=\"clock\" >";
	} else if ($status == 7) {   $clock = "<img src=\"../images/07.png\" alt-text=\"07\" class=\"clock\" >";
	} else if ($status == 8) {   $clock = "<img src=\"../images/08.png\" alt-text=\"08\" class=\"clock\" >";
	} else if ($status == 9) {   $clock = "<img src=\"../images/09.png\" alt-text=\"09\" class=\"clock\" >";
	} else if ($status == 10) {  $clock = "<img src=\"../images/10.png\" alt-text=\"10\" class=\"clock\" >";
	} else if ($status == 11) {  $clock = "<img src=\"../images/11.png\" alt-text=\"11\" class=\"clock\" >";
	} else {                     $clock = "[ * ]";
	}

	$dirFigureBase = $directory."users/".$user."/projects/".$project."/";
	$urlFigureBase = $url."users/".$user."/projects/".$project."/";

	// Load 'dataType' from project folder.
	$handle   = fopen($dirFigureBase."dataType.txt", "r");
	$dataType = trim(fgets($handle));
	fclose($handle);

	// Load 'parent' from project folder.
	$handle = fopen($dirFigureBase."parent.txt", "r");
	$parent = trim(fgets($handle));
	fclose($handle);

	echo "\n<!--\tuser    = ".$user;
	echo "\n\tproject = ".$project." --!>";

	if (file_exists($dirFigureBase."complete.txt")) {
		echo "\n<!-- complete file found.\n--!>";
		if (strcmp($dataType,"0") == 0) {			// SnpCghArray
			$linear_png_url      = $urlFigureBase."fig.CNV-SNP-map.2.png";
			$linear_eps_url      = $urlFigureBase."fig.CNV-SNP-map.2.eps";
			$standard_png_url    = $urlFigureBase."fig.CNV-SNP-map.1.png";
			$standard_eps_url    = $urlFigureBase."fig.CNV-SNP-map.1.eps";
			$CGH_GCbias          = $urlFigureBase."fig_GCratio_vs_CGH.png";
		} else {
			if ((strcmp($dataType[0],"2") == 0) && (strcmp($parent,$project) == 0)) { // ddRADseq, must check if is a parental control dataset => no figure display.
				$isRADseqParent = 1;
			} else {
				$isRADseqParent = 0;
				if (file_exists($dirFigureBase."fig.CNV-LOH-map.2.png")) {
					$linear_png_url      = $urlFigureBase."fig.CNV-LOH-map.2.png";
					$linear_eps_url      = $urlFigureBase."fig.CNV-LOH-map.2.eps";
					$standard_png_url    = $urlFigureBase."fig.CNV-LOH-map.1.png";
					$standard_eps_url    = $urlFigureBase."fig.CNV-LOH-map.1.eps";
				} else if (file_exists($dirFigureBase."fig.CNV-SNP-map.2.png")) {
					$linear_png_url      = $urlFigureBase."fig.CNV-SNP-map.2.png";
					$linear_eps_url      = $urlFigureBase."fig.CNV-SNP-map.2.eps";
					$standard_png_url    = $urlFigureBase."fig.CNV-SNP-map.1.png";
					$standard_eps_url    = $urlFigureBase."fig.CNV-SNP-map.1.eps";
				} else {
					$linear_png_url      = $urlFigureBase."fig.CNV-map.2.png";
					$linear_eps_url      = $urlFigureBase."fig.CNV-map.2.eps";
					$standard_png_url    = $urlFigureBase."fig.CNV-map.1.png";
					$standard_eps_url    = $urlFigureBase."fig.CNV-map.1.eps";
				}
				$bias_figure_dir = $dirFigureBase."fig.examine_bias.png";
				$bias_figure_url = $urlFigureBase."fig.examine_bias.png";
			}
		}
		$CGD_annotations_dir = $dirFigureBase."CGD_annotations.".$project.".txt";
		$CGD_annotations_url = $urlFigureBase."CGD_annotations.".$project.".txt";
		$manualLOH_png_dir   = $dirFigureBase."fig.CNV-manualLOH-map.2.png";
		$manualLOH_png_url   = $urlFigureBase."fig.CNV-manualLOH-map.2.png";
		$manualLOH_eps_url   = $urlFigureBase."fig.CNV-manualLOH-map.2.eps";

		if ($isRADseqParent == 1) {
			?>
			<html>
			<body onload = "parent.resize_iframe('<?php echo $key; ?>', 32);" >
			<div class="tab">
			Dataset is a parental control for CNV normalization during ddRADseq analysis. No figure is constructed.
			</div>
			</body>
			</html>
			<?php
		} else {
			?>
			<html>
			<body onload = "parent.resize_iframe('<?php echo $key; ?>', 115);" >
			<table border="0" cellpadding="5"><tr><td>
			<table border="0" style="max-width:100%" marginwidth="0" marginheight="0" vspace="0" hspace="0"><tr><td valign="middle">
			<div id="imageContainer"><img src="<?php echo $linear_png_url; ?>" style="max-width:100%"></img>
			</td><td align="right" valign="top">
			<font size="2">
			Linear
			<img src="../images/icon_png_15b.png" alt-text="linear, PNG"   onclick="loadImage('<?php echo $linear_png_url;    ?>',100,115);">
			<img src="../images/icon_eps_15b.png" alt-text="linear, EPS"   onclick="loadExternal('<?php echo $linear_eps_url; ?>');"
			<br><br>
			Standard
			<img src="../images/icon_png_15b.png" alt-text="standard, PNG" onclick="loadImage('<?php echo $standard_png_url;    ?>',50 ,440);">
			<img src="../images/icon_eps_15b.png" alt-text="standard, EPS" onclick="loadExternal('<?php echo $standard_eps_url; ?>');"
			<br><br>
			<?php
			if (file_exists($manualLOH_png_dir)) {
				?>
				mLOH
				<img src="../images/icon_png_15b.png" alt-text="manual LOH, linear, PNG"   onclick="loadImage('<?php echo $manualLOH_png_url;    ?>',100,115);">
				<img src="../images/icon_eps_15b.png" alt-text="manual LOH, linear, EPS"   onclick="loadExternal('<?php echo $manualLOH_eps_url; ?>');"
				<br><br>
				<?php
			}
			if (file_exists($CGD_annotations_dir)) {
				echo "<button onclick=\"loadExternal('".$CGD_annotations_url."',50,  220);\">GBrowse</button><br>";
			}
			if (file_exists($bias_figure_dir)) {
				echo "<button onclick=\"loadExternal('".$bias_figure_url."',50,  220);\">read bias</button><br>";
			}
			?>
			</font>
			</td></tr></table>
			</td></tr></table>
			</body>
			</html>
			<?php
		}
	} else if (file_exists($dirFigureBase."error.txt")) {
		echo "\n<!-- error file found.\n--!>";
		// Load error.txt from project folder.
        $handle = fopen($dirFigureBase."error.txt", "r");
        $error = fgets($handle);
        fclose($handle);
		?>
		<html>
		<body onload = "parent.resize_iframe('<?php echo $key; ?>', 115);" >
			<font color="red"><b>[Error : Consult site admin.]</b></font><br>
			<?php echo $error; ?>
		</body>
		</html>
		<?php
	} else if (file_exists($dirFigureBase."working.txt")) {
		echo "\n<!-- working file found. --!>\n";
		// Load last line from "condensed_log.txt" file.
		$condensedLog      = explode("\n", trim(file_get_contents($dirFigureBase."condensed_log.txt")));
		$condensedLogEntry = $condensedLog[count($condensedLog)-1];
		?>
		<script type="text/javascript">
		var user    = "<?php echo $user; ?>";
		var project = "<?php echo $project; ?>";
		var key     = "<?php echo $key+1; ?>";
		var status  = "<?php echo $status; ?>";
		reload_page=function() {
			// Make a form to generate a form to POST information to pass along to page reloads, auto-triggered by form submit.
			var autoSubmitForm = document.createElement('form');
			    autoSubmitForm.setAttribute('method','post');
			    autoSubmitForm.setAttribute('action','project.working_server.php');
			var input2 = document.createElement('input');
			    input2.setAttribute('type','hidden');
			    input2.setAttribute('name','key');
			    input2.setAttribute('value',key);
			    autoSubmitForm.appendChild(input2);
			var input2 = document.createElement('input');
			    input2.setAttribute('type','hidden');
			    input2.setAttribute('name','user');
			    input2.setAttribute('value',user);
			    autoSubmitForm.appendChild(input2);
			var input3 = document.createElement('input');
			    input3.setAttribute('type','hidden');
			    input3.setAttribute('name','project');
			    input3.setAttribute('value',project);
			    autoSubmitForm.appendChild(input3);
			var input4 = document.createElement('input');
			    input4.setAttribute('type','hidden');
			    input4.setAttribute('name','status');
			    input4.setAttribute('value',status);
			    autoSubmitForm.appendChild(input4);
			autoSubmitForm.submit();
		}
		// Initiate recurrent call to reload_page function, which depends upon project status.
		var internalIntervalID = window.setInterval(reload_page, 3000);
		</script>
		<html>
		<body onload = "parent.resize_iframe('<?php echo $key; ?>', 20*2+12);" class="tab">
		<font color="red"><b>[Processing uploaded data.]</b></font>
		<?php
		echo $clock."<br>";
		if (strcmp($dataType,"0") == 0) {
			echo "SnpCgh microarray analysis usually complete in a few minutes.";
		} else {
			echo $condensedLogEntry;
		}
	} else {
		echo "\n<html>\n<body>\n";
		echo "dirBase = ".$dirFigureBase."<br>\n";
		echo "urlBase = ".$urlFigureBase."<br>\n";
		echo "complete.txt file not found properly.<br>\n";
	}
	echo "\n";
	?>
	</body>
	<script type="text/javascript">
		function loadImage(imageUrl,imageScale,iframeHeight) {
			document.getElementById('imageContainer').innerHTML = "<img src=\""+imageUrl+"\" style=\"max-width:"+imageScale+"%\"></img>";
			parent.resize_iframe('<?php echo $key; ?>', iframeHeight);
		}
		function loadExternal(imageUrl) {
			window.open(imageUrl);
		}
	</script>
</HTML>
