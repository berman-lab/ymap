<?php

// return the current size in GB of the user folder
function getUserUsageSize($userName)
{
	// get all finished folders concatenated in one line in order to call du
	$finishedFolders = shell_exec("find " . "users/".$userName . "/  -type f -iname 'complete.txt' | sed -e \"s/complete.txt//g\" | tr '\n' ' ' ");
	// check if finished folders were found if so calculating the size else return 0
	if ($finishedFolders != "")
		return shell_exec("du -scm " . $finishedFolders ."| awk 'END{print $1}'") / (1000);
	return 0;
}

// return the size of the user quota in GB
function getUserQuota($userName)
{
	// load hardcoded quota from constants
	require_once 'constants.php';
	// check if user has a personal quota if so overriding quota
	if (file_exists("users/".$userName . "/quota.txt"))
		$quota = trim(file_get_contents("users/".$userName . "/quota.txt"));
	return $quota;
}

?>