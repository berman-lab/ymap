<?php
	session_start();
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<HTML>
<HEAD>
	<style type="text/css">
		.upload {
			width:          675px;      // 675px;
			border:         0;
			height:         40px;   // 40px;
			vertical-align: middle;
			align:          left;
			margin:         0px;
			overflow:       hidden;
		}
		html, body {
			margin:         0px;
			border:         0;
			overflow:       hidden;
		}
	</style>
<meta http-equiv="content-type" content="text/html; charset=iso-8859-1">
<title>Install project into pipeline.</title>
</HEAD>
<?php
    require_once 'constants.php';
	$fileName        = filter_input(INPUT_POST, "fileName",        FILTER_SANITIZE_STRING);
	$user            = filter_input(INPUT_POST, "user",            FILTER_SANITIZE_STRING);
	$project         = filter_input(INPUT_POST, "project",         FILTER_SANITIZE_STRING);
	$key             = filter_input(INPUT_POST, "key",             FILTER_SANITIZE_STRING);

// Generate 'datafile1.txt' file containing: name of first data file.
	$outputName = $directory."users/".$user."/projects/".$project."/datafile1.txt";
	$output     = fopen($outputName, 'w');
	fwrite($output, $fileName);
	fclose($output);
	unset($outputName);
	unset($output);

?>
<BODY onload = "parent.resize_iframe('<?php echo $key; ?>', 80);" >
<div id="project_frame"></div>

<script type="text/javascript">
onload=function(){
	// Load 'uploader.html' into an internal frame for each project needing a file uploaded.
	var p_el = document.getElementById("project_frame");
	p_el.innerHTML="<iframe id=\"<?php echo $key ?>\" name=\"<?php echo $key ?>\" class=\"upload\" src=\"<?php echo $url; ?>uploader.html\">";
	// Load file destination and concluding processing script into internal frame.   ddd
	top.frames['<?php echo $key ?>'].frames[0].label             = "Exp FASTQ (2/2)...";
	top.frames['<?php echo $key ?>'].frames[0].target_dir        = "<?php echo $directory."users/".$user."/projects/".$project."/"; ?>";
	top.frames['<?php echo $key ?>'].frames[0].conclusion_script = "<?php echo $url."php/ddRADseq_4.install_2.php";                 ?>";
	top.frames['<?php echo $key ?>'].frames[0].user              = "<?php echo $user;                                               ?>";
	top.frames['<?php echo $key ?>'].frames[0].project           = "<?php echo $project;                                            ?>";
	top.frames['<?php echo $key ?>'].frames[0].key               = "<?php echo $key;                                                ?>";
}
</script>

</BODY>
</HTML>
