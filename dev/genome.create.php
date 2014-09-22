<?php
	session_start();
	require_once 'php/constants.php';
	// If the user is not logged on, redirect to login page.
//	if(!isset($_SESSION['logged_on'])){
//		header('Location: user.login.php');
//	} else {
		echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";
//	}
?>
<html lang="en">
	<head>
		<style type="text/css">
			body {font-family: arial;}
		</style>
		<meta http-equiv="content-type" content="text/html; charset=utf-8">
		<title>[Needs Title]</title>
	</head>
	<body>
		<div id="loginControls"><p>
			<?php if(isset($_SESSION['logged_on'])){echo "user: ".$_SESSION['user']."<br>";}?>

			<!-- If user is logged in, show logout button, otherwise, show the login button so we can get the user logged in-->
			<button type="button" onclick="window.location.href='<?php if(isset($_SESSION['logged_on'])){echo "php/logout_server.php";}else{echo "user.login.php";}?>'"><?php if(isset($_SESSION['logged_on'])){echo "Logout";}else{echo "Login";}?></button>
			<button type="button" onclick="window.location.href='index.php'">Back to Home</button>
		</p></div>
		<div id="genomeCreationInformation"><p>
			<form action="php/genome.create_server.php" method="post">
				<label for="newGenomeName">Genome Name: </label><input type="text" name="newGenomeName" id="newGenomeName" size="100"><br>
				<br>
				<input type="submit" value="Create New Genome">
			</form>
		</p></div>
	</body>
</html>
