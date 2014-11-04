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

// Generate 'working.txt' file to let pipeline know processing is started.
	$outputName = $directory."users/".$user."/projects/".$project."/working.txt";
	$output     = fopen($outputName, 'w');
	fwrite($output, "working");
	fclose($output);
	chmod($outputName,0644);

// Generate 'datafiles.txt' file containing: name of first data file.
	$outputName = $directory."users/".$user."/projects/".$project."/datafiles.txt";
	$output     = fopen($outputName, 'w');

	$fileNames = explode(",", $fileName);
	foreach ($fileNames as $name) {
		$name  = str_replace("\\", ",", $name);
		fwrite($output, $name."\n");
	}

	fclose($output);
	unset($outputName);
	unset($output);

?>
<BODY onload = "parent.resize_iframe('<?php echo $key; ?>', 40);" >
<font color="red">[Four file upload completed.]<br>[Processing...]</font>
<?php
	echo "<br>";
	print_r($fileNames);
?>
</BODY>
</HTML>

<script type="text/javascript">
//<BODY>
//<font size="2" color="red">Upload complete; processing...</font><br>
//<script type="text/javascript">
//	var user    = "<?php echo $user;    ?>";
//	var project = "<?php echo $project; ?>";
//	var key     = "<?php echo $key;     ?>";
//	// construct and submit form to move on to "project.working_server.php";
//	var autoSubmitForm = document.createElement("form");
//		autoSubmitForm.setAttribute("method","post");
//		autoSubmitForm.setAttribute("action","project.working_server.php");
//	var input2 = document.createElement("input");
//		input2.setAttribute("type","hidden");
//		input2.setAttribute("name","key");
//		input2.setAttribute("value",key);
//		autoSubmitForm.appendChild(input2);
//	var input2 = document.createElement("input");
//		input2.setAttribute("type","hidden");
//		input2.setAttribute("name","user");
//		input2.setAttribute("value",user);
//		autoSubmitForm.appendChild(input2);
//	var input3 = document.createElement("input");
//		input3.setAttribute("type","hidden");
//		input3.setAttribute("name","project");
//		input3.setAttribute("value",project);
//		autoSubmitForm.appendChild(input3);
//	autoSubmitForm.submit();
//</script>
//</BODY>
//</HTML>
</script>
