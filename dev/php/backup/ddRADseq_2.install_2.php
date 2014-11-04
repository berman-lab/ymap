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
	$fileName = filter_input(INPUT_POST, "fileName", FILTER_SANITIZE_STRING);
	$user     = filter_input(INPUT_POST, "user",     FILTER_SANITIZE_STRING);
	$project  = filter_input(INPUT_POST, "project",  FILTER_SANITIZE_STRING);
	$key      = filter_input(INPUT_POST, "key",      FILTER_SANITIZE_STRING);

// Generate 'datafile2.txt' file containing: name of second data file.
	$outputName = $directory."users/".$user."/projects/".$project."/datafile2.txt";
	$output     = fopen($outputName, 'w');
	fwrite($output, $fileName);
	fclose($output);
	unset($outputName);
	unset($output);
?>
<BODY onload = "parent.resize_iframe('<?php echo $key; ?>', 40);" >
<font color="red">[Two file installation completed.]<br>[Processing...]</font>
</BODY>
</HTML>
