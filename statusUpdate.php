<?php
	session_start();
	require_once 'constants.php';
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
	// Load last line from "status_log.txt" file.
	$statusLogLines = explode("\n", trim(file_get_contents("status_log.txt")));
	$statusLog      = str_replace("\n","<br>\n",trim(file_get_contents("status_log.txt")));

	if (count($statusLogLines) > 1) {
		$statusLogEntry = $statusLog;
	} else if (count($statusLogLines) == 1) {
		if ($statusLog[count($statusLog)-1] == "") {
			$statusLogEntry = "No system status notes.";
		} else {
			$statusLogEntry = $statusLog;
		}
	} else {
		$statusLogEntry = "No system status notes.";
	}
?>
<script type="text/javascript">
reload_page=function() {
	var autoSubmitForm = document.createElement('form');
	    autoSubmitForm.setAttribute('method','post');
	    autoSubmitForm.setAttribute('action','statusUpdate.php');
	    autoSubmitForm.submit();
}
var internalIntervalID = window.setInterval(reload_page, 30000);
</script>
<body onload = "parent.resize_iframe('status', 75);" class="tab">
	<font size="2">
	<?php echo $statusLogEntry; ?>
	</font>
</BODY>
</HTML>
