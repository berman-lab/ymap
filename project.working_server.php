<?php
	session_start();
	$user = $_SESSION['user'];
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

	$bad_chars = array("~","@","#","$","%","^","&","*","(",")","+","=","|","{","}","<",">","?",".",",","\\","/","'",'"',"[","]","!");
	$project   = str_replace($bad_chars,"",trim(filter_input(INPUT_POST, "project", FILTER_SANITIZE_STRING)));
	$key       = str_replace($bad_chars,"",trim(filter_input(INPUT_POST, "key", FILTER_SANITIZE_STRING)));
	$status   = filter_input(INPUT_POST, "status",   FILTER_SANITIZE_STRING);

	// increment clock animation...
	$status   = ($status + 1) % 12;
	if ($status == 0) {          $clock = "<img src=\"images/12.png\" alt-text=\"12\" class=\"clock\" >";
	} else if ($status == 1) {   $clock = "<img src=\"images/01.png\" alt-text=\"01\" class=\"clock\" >";
	} else if ($status == 2) {   $clock = "<img src=\"images/02.png\" alt-text=\"02\" class=\"clock\" >";
	} else if ($status == 3) {   $clock = "<img src=\"images/03.png\" alt-text=\"03\" class=\"clock\" >";
	} else if ($status == 4) {   $clock = "<img src=\"images/04.png\" alt-text=\"04\" class=\"clock\" >";
	} else if ($status == 5) {   $clock = "<img src=\"images/05.png\" alt-text=\"05\" class=\"clock\" >";
	} else if ($status == 6) {   $clock = "<img src=\"images/06.png\" alt-text=\"06\" class=\"clock\" >";
	} else if ($status == 7) {   $clock = "<img src=\"images/07.png\" alt-text=\"07\" class=\"clock\" >";
	} else if ($status == 8) {   $clock = "<img src=\"images/08.png\" alt-text=\"08\" class=\"clock\" >";
	} else if ($status == 9) {   $clock = "<img src=\"images/09.png\" alt-text=\"09\" class=\"clock\" >";
	} else if ($status == 10) {  $clock = "<img src=\"images/10.png\" alt-text=\"10\" class=\"clock\" >";
	} else if ($status == 11) {  $clock = "<img src=\"images/11.png\" alt-text=\"11\" class=\"clock\" >";
	} else {                     $clock = "[ * ]";
	}

	$dirFigureBase = "users/".$user."/projects/".$project."/";

	// Load 'dataFormat' from project folder.
	$handle   = fopen($dirFigureBase."dataFormat.txt", "r");
	$dataFormat = trim(fgets($handle));
	fclose($handle);

	// Load 'parent' from project folder.
	$handle = fopen($dirFigureBase."parent.txt", "r");
	$parent = trim(fgets($handle));
	fclose($handle);

	echo "\n<!--\tuser    = ".$user;
	echo "\n\tproject = ".$project." --!>";

	$sizeString_1 = "";
	$sizeString_2 = "";

	$sizeFile_1   = "users/".$user."/projects/".$project."/upload_size_1.txt";
	$handle       = @fopen($sizeFile_1,'r');
	if ($handle) {
		$sizeString_1 = trim(fgets($handle));
		fclose($handle);
	}

	$sizeFile_2   = "users/".$user."/projects/".$project."/upload_size_2.txt";
	$handle       = @fopen($sizeFile_2,'r');
	if ($handle) {
		$sizeString_2 = trim(fgets($handle));
		fclose($handle);
	}

	if (($sizeString_1 !== "") || ($sizeString_2 !== "")) {
		echo "\n<script type='text/javascript'>\n";
		echo "parent.parent.update_project_file_size('".$key."','".$sizeString_1."','".$sizeString_2."');";
		echo "\n</script>\n";
	}

	if (file_exists($dirFigureBase."complete.txt")) {
		echo "\n<!-- complete file found.\n--!>";
		// Hide iframe and adjust color of entry to indicate completion.
		?>
		<html>
		<body onload = "parent.parent.update_project_label_color('<?php echo $key; ?>','#00AA00'); parent.parent.resize_project('<?php echo $key; ?>', 0); parent.parent.update_project_remove_iframe('<?php echo $key; ?>', '<?php echo htmlspecialchars(json_encode(scandir("users/$user/projects/$project"))); ?>');">
		</body>
		</html>
		<?php
	} else if (file_exists($dirFigureBase."error.txt")) {
		echo "\n<!-- error file found.\n--!>";
		// Load error.txt from project folder.
        $handle = fopen($dirFigureBase."error.txt", "r");
        $error = fgets($handle);
        fclose($handle);
		?>
		<html>
		<body onload = "parent.parent.resize_project('<?php echo $key; ?>', 100);" >
			<font color="red"><b>[Error : Consult site admin.]</b></font><br>
			<?php echo $error; ?>
		</body>
		</html>
		<?php
	} else if (file_exists($dirFigureBase."working.txt")) {
		// Load start time from 'working.txt'
		$startTimeStamp = file_get_contents("users/".$user."/projects/".$project."/working.txt");
		$startTime      = strtotime($startTimeStamp);
		$currentTime    = time();
		$intervalTime   = $currentTime - $startTime;
		if (strcmp($dataFormat,"0") == 0) {
			$timeLimit  = 60*60;   // 1 hour minutes for array analysis.
		} else {
			$timeLimit  = 60*60*6; // 6 hours
		}
		if ($intervalTime > $timeLimit) { // likely error.
			?>
			<BODY onload = "parent.parent.resize_project('<?PHP echo $key; ?>', 100);" class="tab">
				<font color="red" size="2"><b>[Error]</b></font><?php echo " &nbsp; &nbsp; ".$clock; ?><br>
				<font size="2">Processing of data has taken longer than expected and might be stalled.<br>
				Check back later or contact the admin through the "System" tab with details and they will check on the job.<br>
				Don't delete the job until the admin has responded, or they will be unable to assist.</font>
			</BODY>
			</HTML>
			<?PHP
		} else {
			echo "\n<!-- working file found. --!>\n";
			// Load last line from "condensed_log.txt" file.
			$condensedLog      = explode("\n", trim(file_get_contents($dirFigureBase."condensed_log.txt")));
			$condensedLogEntry = $condensedLog[count($condensedLog)-1];
			?>
			<script type="text/javascript">
			var user    = "<?php echo $user; ?>";
			var project = "<?php echo $project; ?>";
			var key     = "<?php echo $key; ?>";
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
				document.body.appendChild(autoSubmitForm);
				autoSubmitForm.submit();
			}
			// Initiate recurrent call to reload_page function, which depends upon project status.
			var internalIntervalID = window.setInterval(reload_page, 3000);
			</script>
			<html>
			<body onload = "parent.parent.resize_project('<?PHP echo $key; ?>', 38); parent.parent.update_project_label_color('<?php echo $key; ?>','#BB9900');" class="tab">
			<div style='color: red; display: inline-block; font-size: 10px;'><b>[Processing uploaded data.]</b></div>
			<?php
			echo $clock."<br>";
			if (strcmp($dataFormat,"0") == 0) {
				echo "<div style='font-size: 10px;'>";
				echo "SnpCgh microarray analysis usually completes in less than an hour, depending on system load.";
				echo "</div>";
			} else {
				echo "<div style='font-size: 10px;'>";
				echo $condensedLogEntry;
				echo "</div>";
			}
		}
	} else {
		echo "\n<html>\n<body>\n";
		echo "complete.txt file not found properly.<br>\n";
	}
	echo "\n";
	?>
	</body>
</HTML>
