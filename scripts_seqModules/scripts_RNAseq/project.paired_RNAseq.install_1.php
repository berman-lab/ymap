<?php
	session_start();
	error_reporting(E_ALL);
        require_once '../../constants.php';
	require_once '../../POST_validation.php';
        ini_set('display_errors', 1);

        // If the user is not logged on, redirect to login page.
        if(!isset($_SESSION['logged_on'])){
		session_destroy();
                header('Location: ../../');
        }

	// load strings from session.
	$user     = $_SESSION['user'];
	$fileName = $_SESSION['fileName'];
	$project  = $_SESSION['project'];
	$key      = $_SESSION['key'];

	$project_dir = "../../users/".$user."/projects/".$project;
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<HTML>
<HEAD>
	<style type="text/css">
		body {font-family: arial;}
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
// Initialize log file.
	$logOutputName = $project_dir."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'w');
	fwrite($logOutput, "Log file initialized.\n");
	fwrite($logOutput, "#..............................................................................\n");
	fwrite($logOutput, "Running 'php/project.paired_RNAseq.install_1.php'.\n");
	fwrite($logOutput, "Variables passed via POST :\n");
	fwrite($logOutput, "\tfileName = '".$fileName."'\n");
	fwrite($logOutput, "\tuser     = '".$user."'\n");
	fwrite($logOutput, "\tproject  = '".$project."'\n");
	fwrite($logOutput, "\tkey      = '".$key."'\n");
	fwrite($logOutput, "#============================================================================== 1\n");

	$condensedLogOutputName = $project_dir."/condensed_log.txt";
	$condensedLogOutput     = fopen($condensedLogOutputName, 'w');
	fwrite($condensedLogOutput, "Initializing.\n");
	fclose($condensedLogOutput);
	chmod($outputName,0755);

// Generate 'working.txt' file to let pipeline know processing is started.
	$outputName      = $project_dir."/working.txt";
	$output          = fopen($outputName, 'w');
	$startTimeString = date("Y-m-d H:i:s");
	fwrite($output, $startTimeString);
	fclose($output);
	chmod($outputName,0755);
	fwrite($logOutput, "\tGenerated 'working.txt' file.\n");

// Installation continues with next php script... strings recieved as POST are forwarded to next script.
	fwrite($logOutput, "Passing control to : 'php/project.paired_RNAseq.install_2.php'\n");
 	$system_call_string = "php project.paired_RNAseq.install_2.php ".$fileName." ".$user." ".$project." > /dev/null &";
	system($system_call_string);

	fwrite($logOutput, "Loading 'project.working_server.php' into iframe.\n");
	fclose($logOutput);
?>
<font size="2" color="red">Upload complete; processing...</font><br>
<script type="text/javascript">
	var user    = "<?php echo $user;    ?>";
	var project = "<?php echo $project; ?>";
	var key     = "<?php echo $key;     ?>";
	var status  = "0";
	// construct and submit form to move on to "project.working_server.php";
	var autoSubmitForm = document.createElement("form");
		autoSubmitForm.setAttribute("method","post");
		autoSubmitForm.setAttribute("action","../project.working_server.php");
	var input2 = document.createElement("input");
		input2.setAttribute("type","hidden");
		input2.setAttribute("name","key");
		input2.setAttribute("value",key);
		autoSubmitForm.appendChild(input2);
	var input2 = document.createElement("input");
		input2.setAttribute("type","hidden");
		input2.setAttribute("name","user");
		input2.setAttribute("value",user);
		autoSubmitForm.appendChild(input2);
	var input3 = document.createElement("input");
		input3.setAttribute("type","hidden");
		input3.setAttribute("name","project");
		input3.setAttribute("value",project);
		autoSubmitForm.appendChild(input3);
	var input4 = document.createElement("input");
		input4.setAttribute("type","hidden");
		input4.setAttribute("name","status");
		input4.setAttribute("value",status);
	document.body.appendChild(autoSubmitForm);
	autoSubmitForm.submit();
</script>
</BODY>
</HTML>
