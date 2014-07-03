function deleteProjectConfirmation(project,key){
	$("#p2_"+project+"_delete").html("<b><font color=\"red\">[Are you sure?]</font><button type='button' onclick='deleteProject_yes(\""+project+"\",\""+key+"\")'>Yes, delete.</button><button type='button' onclick='deleteProject_no(\""+project+"\",\""+key+"\")'>No, cancel</button></b>");
	document.getElementById(key).style.display = 'none';
}

// Requires javascript value user to be defined globally
function deleteProject_yes(project,key){
	$.ajax({
		url : 'php/project.delete_server.php',
		type : 'post',
		data : {
			user: user,
			project: project
		},
		success : function(answer){
			if(answer == "COMPLETE"){
				window.location.href=window.location.href;
			}
		}
	});
	update_projectsShown_after_project_delete(key);
}
function deleteProject_no(project,key){
	$("#p2_"+project+"_delete").html("");
	document.getElementById(key).style.display = 'inline';
}

