%Updates the textbox that displays the Synapse Mode
function updatemode(app) 
    switch app.syn.getMode
        case 0
            app.CurrentModeLabel.Text = 'Idle';
        case 1
            app.CurrentModeLabel.Text = 'Standby';
        case 2
            app.CurrentModeLabel.Text = 'Preview';
        case 3
            app.CurrentModeLabel.Text = 'Recording';
    end
end