<?php
	session_start();
	$user = $_SESSION['user'];

	$project_key       = filter_input(INPUT_GET, "key",      FILTER_SANITIZE_STRING);
	$user              = filter_input(INPUT_GET, "user",     FILTER_SANITIZE_STRING);
	$project           = filter_input(INPUT_GET, "project",  FILTER_SANITIZE_STRING);
	$conclusion_script = filter_input(INPUT_GET, "script",   FILTER_SANITIZE_STRING);
?>
<!DOCTYPE html>
<head>
<body>
<?php
$FTPdrop_dir   = "FTP_drop/".$user."/";
$FTPdrop_files = array_diff(glob($FTPdrop_dir."*"), array('..', '.'));
// Sort files by date, newest first.
array_multisort(SORT_DESC, $FTPdrop_files);
// Trim path from each folder string.
foreach($FTPdrop_files as $key=>$file) {
	$FTPdrop_files[$key] = str_replace($FTPdrop_dir,"",$file);
}
$FTPdrop_fileCount = count($FTPdrop_files);
?>
<form action="file_selection.2.server.php" method="post">
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Select datafile <select name="File_selection_1">
<?php
foreach($FTPdrop_files as $key=>$file) {
	echo "<option value=\"".$key."\">".$FTPdrop_files[$key]."</option>\n";
}
?>
</select> and <select name="File_selection_2">
<?php
foreach($FTPdrop_files as $key=>$file) {
	echo "<option value=\"".$key."\">".$FTPdrop_files[$key]."</option>\n";
}
?>
</select>
to
<input type="submit" id="submitbutton" value="load into project" onclick="document.getElementById('submitbutton').style.display='none'; document.getElementById('submitbutton_replace').style.display='inline';">
<div id="submitbutton_replace" style="display:none">load into project</div>.
<input type="hidden" name="key"        value="<?php echo $project_key;       ?>">
<input type="hidden" name="user"       value="<?php echo $user;              ?>">
<input type="hidden" name="project"    value="<?php echo $project;           ?>">
<input type="hidden" name="script"     value="<?php echo $conclusion_script; ?>">
</form>

<p id="demo" onclick="myFunction()">Click me to change my text color.</p>
<script>
function myFunction() {
    document.getElementById("demo").style.color = "red";
}
</script>

</body>
</html>
