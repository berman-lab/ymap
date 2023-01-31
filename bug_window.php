<?php
	session_start();
	if(!isset($_SESSION['logged_on'])){ ?> <script type="text/javascript"> parent.reload(); </script> <?php } else { $user = $_SESSION['user']; }
	require_once 'constants.php';
	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";
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
					description = $('#description').val();
					$.ajax({
						url : 'addBug_server.php',
						type : 'post',
						data : {
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
			Describe your bug or feature request below.<br>
			If you are having a problem with a specific project, make sure to include the project name in your comment.<br>
			If you are responding to a previous comment, make sure to include the comment number from the first collumn.<br>
			Once an admin has responded, the comment entries will be colored <font style='background-color: #AAFFAA'>green to
			indicate a resolved issue</font>, <font style='background-color: #FFFFAA'>yellow to indicate an ongoing issue</font>,
			or <font style='background-color: #FFAAAA'>red to indicate the need for more information/discussion</font>.<br>
			Your comments will only be visible to you and site administrators.<br>
			<textarea id="description" rows="5" cols="150">Comment here.</textarea><br />
			<button type="button" onclick="submitBug()">Submit Bug or Feature!</button><br /><br />
		</div>
		<div id="bugList">
			<table id="bugTable" border="1">
				<tr bgcolor="#DDDDDD">
					<th>#</th>
					<th>User</th>
					<th colspan=80>Description</th>
					<th>Type</th>
					<th colspan=80>Admin Notes</th>
				</tr>
					<?php
					$bugtrackFile = "bugtracker/bugs.txt";
					if (file_exists($bugtrackFile)) {
						$bugFileHandle = fopen($bugtrackFile, 'r');
						$comment_count = 0;
						while (($buffer = fgets($bugFileHandle, 4096)) !== false) {
							if (strcmp($buffer,"") <> 1) {   $comment_count = $comment_count + 1;   }
							list ($user, $description, $type, $notes, $status) = explode("|", $buffer);
							$super_user_flag_file = "users/".$user."/super.txt";
							if (file_exists($super_user_flag_file)) {  // Super-user privilidges.
								if ($status == 1) {              echo "<tr bgcolor='FFAAAA'>";
								} else if ($status == 2) {       echo "<tr bgcolor='FFFFAA'>";
								} else if ($status == 3) {       echo "<tr bgcolor='AAFFAA'>";
								} else {                         echo "<tr>";
								}
								if (strcmp($buffer,"") == 1) {   echo "<td align='center'></td>";
								} else {                         echo "<td align='center'>".$comment_count."</td>";
								}
								echo "<td align='center'>".$user."</td>\n<td colspan=80>".$description."</td>\n<td align='center'>".$type."</td>\n<td colspan=80>".$notes."</td>\n";
								echo "</tr>";
							} else {
								if ($user == 'admin') {
									if ($status == 1) {              echo "<tr bgcolor='FFAAAA'>";
									} else if ($status == 2) {       echo "<tr bgcolor='FFFFAA'>";
									} else if ($status == 3) {       echo "<tr bgcolor='AAFFAA'>";
									} else {                         echo "<tr>";
									}
									if (strcmp($buffer,"") == 1) {   echo "<td align='center'></td>";
									} else {                         echo "<td align='center'>".$comment_count."</td>";
									}
									echo "<td align='center'>".$user."</td>\n<td colspan=80>".$description."</td>\n<td align='center'>".$type."</td>\n<td colspan=80>".$notes."</td>\n";
									echo "</tr>";
								}
							}
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
