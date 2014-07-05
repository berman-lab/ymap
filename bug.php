<?php
    session_start();
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
        "http://www.w3.org/TR/html4/loose.dtd">
<?php
	require_once 'php/constants.php';
	chmod($GLOBALS['directory'].$GLOBALS['bugtrackDir'].$GLOBALS['bugtrackFile'], 0777);

	$user = $_SESSION['user'];
?>
<html lang="en">
	<head>
		<style type="text/css">
			body {font-family: arial;}
		</style>
		<meta http-equiv="content-type" content="text/html; charset=utf-8">
		<title>Ymap | Bug Report or Feature Request</title>
		<script type="text/javascript" src="js/jquery-1.7.2.js"></script>
		<script type="text/javascript" src="js/jquery.form.js"></script>
		<script>
				function submitBug(project){
					submitter = $('#submitter').val();
					description = $('#description').val();
					$.ajax({
						url : '<?php echo $GLOBALS['url']?>php/addBug_server.php',
						type : 'post',
						data : {
							submitter: submitter,
							description: description
						},
						success : function(answer){
							window.location.href=window.location.href;
						}
					});
				}
		</script>
	</head>
	<body>
		<div id="bugSubmit">
			<button type="button" onclick="window.location.href='index.php'">Back to Home</button><b>User: <?php echo $user?></b>
			<input type="hidden" id="submitter" name="submitter" value="<?php echo $user?>"><br>
			<textarea id="description" rows="15" cols="150">Describe your bug or feature request. If you are having a problem with a specific project, make sure to include the project name in your comment. Your comments will only be visible to you and site administrators.</textarea><br />
			<button type="button" onclick="submitBug()">Submit Bug or Feature!</button><br /><br />
		</div>
		<div id="bugList">
			<table id="bugTable" border="1">
				<tr bgcolor="#DDDDDD">
					<th>Submitter</th>
					<th colspan=80>Description</th>
					<th>Type</th>
					<th colspan=80>Admin Notes</th>
				</tr>
					<?php
					$bugtrackFile = $GLOBALS['directory'].$GLOBALS['bugtrackDir'].$GLOBALS['bugtrackFile'];
					if (file_exists($bugtrackFile)) {
						$bugFileHandle = fopen($bugtrackFile, 'r');
						while (($buffer = fgets($bugFileHandle, 4096)) !== false) {
							list ($submitter, $description, $type, $notes, $status) = explode("|", $buffer);
							$super_user_flag_file = $GLOBALS['directory']."users/".$user."/super.txt";
							if ($status == 1) {
								echo "<tr bgcolor='FFDDDD'>";
							} else if ($status == 2) {
								echo "<tr bgcolor='DDFFDD'>";
							} else {
								echo "<tr>";
							}
							if (file_exists($super_user_flag_file)) {  // Super-user privilidges.
								echo "<td align='center'>".$submitter."</td>
									  <td colspan=80>".$description."</td>
									  <td align='center'>".$type."</td>
									  <td colspan=80>".$notes."</td>";
							} else {
								if (($submitter == $user) || ($submitter == 'admin')) {
									echo "<td align='center'>".$submitter."</td>
										  <td colspan=80>".$description."</td>
										  <td align='center'>".$type."</td>
										  <td colspan=80>".$notes."</td>";
								}
							}
							echo "</tr>";
						}
						fclose($bugFileHandle);
					} else {
						$bugsOutput = fopen($bugtrackFile, 'w');
						chmod($bugsOutput, 0777);
						fwrite($bugtrackFile, "test|test|test|test|0");
						fclose($bugsOutput);
					}
					?>
			</table>
		</div>
	</body>
</html>
