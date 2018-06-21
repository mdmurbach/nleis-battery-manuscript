function start_email()

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.starttls.enable','true');

setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','E_mail','...@....com');
setpref('Internet','SMTP_Username', '...@....com');
setpref('Internet','SMTP_Password', 'password');

end