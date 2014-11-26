<?php
	session_start();
	if (isset($_SESSION['logged_on'])) { $user = $_SESSION['user']; } else { $user = 'default'; }
?>
<style type="text/css">
	html * {
		font-family: arial !important;
	}
</style>
<font size='3'>Report bugs and otherwise communicate with the admin.</font><br><br>
<b>Interact with an admin:</b><br>
<div class="tab">
	<?php
	if (isset($_SESSION['user']) != 0) {
		echo "Report a bug! Request a feature!<br>";
		echo "Or just as the admin a question about how to use a certain feature.<br><br>";
		echo "<b><font color='red'>Site administrator has been occupied with thesis writing and related issues.<br>Response times may vary, but are expected to improve.</font></b><br>";
		echo '<input name="button_BigTracker" type="button" value="Click Here..." onclick="parent.show_hidden(\'Hidden_BugTracker\')">';
	} else {
		echo "Log in using the 'User' tab to gain access to the admin.";
	}
	?>
</div><br>

<b>System status:</b><br><div id="frameContainer.status"></div>
<script type="text/javascript">
var el_status = document.getElementById("frameContainer.status");
<?php
// Load last line from "status_log.txt" file.
$statusLogLines = explode("\n", trim(file_get_contents("status_log.txt")));
$statusLog      = str_replace("\n","<br>",trim(file_get_contents("status_log.txt")));
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
el_status.innerHTML = '<div class="tab"><font size="2"><?php echo $statusLogEntry; ?></font></div>';
</script>
