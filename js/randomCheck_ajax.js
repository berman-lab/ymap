function randomCheck(project, file){
	$.ajax({
		url : 'randomCheck_server.php',
		type : 'post',
		data : {},
		success : function(answer){
			errors = eval(answer);
			
			str = "";
			for(i in errors){
				str = str + errors[i] + "\n";
			}
			str = str + "Alright, so, I just passed the project and file into this function then ran a perl script that outputs random funny stuff. \n";
			str = str + "This can be used to do file validation checks before we allow people to run the script or it could allow us to make files run the scripts more faster by running a sort on them or whatever.\n";
			str = str + "Also, while I used perl anything that can be run and passed arguments on the command line can be used.\n\n";
			
			str = str + "While errors are random at this point, this can be refined once we get the files up there.\n"
			str = str + "Also, having them in an alert is bad since chrome doesn't have a scrollbar on alerts. Later Fix\n";
			str = str + "It is possible to run this immediately after uploading files, that is a thing that can happen.\n";
			
			alert(str);
		}
	});
}
