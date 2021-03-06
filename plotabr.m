function plotabr(app)

ax = [app.Cycle1Axes, app.Cycle2Axes];
y = [app.data.meansig1, app.data.meansig2];
t = (1:size(y, 1)).*1000./app.data.fs;
cololist = [0.72,0.27,1.00; 1, 0, 0; 1.00,0.41,0.16; 0.39,0.83,0.07; 0.00,0.80,0.80; 0.00,0.00,1.00];
tsrow = [str2double(app.Cycle1TimestampRadio.SelectedObject.Text), str2double(app.Cycle2TimestampRadio.SelectedObject.Text)];

for i = 1:2
    tablename = ['TS' num2str(i) 'Table'];

    tstable = app.(tablename).Data;
    for ts = 1:12
        if ~isnan(tstable(ts, 1))
            tstable(ts, 2) = y(round(tstable(ts, 1)/1000*app.data.fs), i);
        end
    end
    app.(tablename).Data = tstable;

    plot(ax(i), 0, 0);
    hold(ax(i), 'on');
    if app.SingleCyclesCheckBox.Value
        plot(ax(i), t, sttrace, 'Color', [0.8 0.8 0.8]);
    end
    plot(ax(i), t, y(:,i), 'k');
    for ts = 1:12
        if ts == tsrow(i)
            plot(ax(i), app.(tablename).Data(ts,1), app.(tablename).Data(ts,2), 'o', 'Color', cololist(rem(ts,6)+1,:),'LineWidth', 2)
            text(ax(i), app.(tablename).Data(ts,1)+.1, app.(tablename).Data(ts,2), num2str(ts), 'FontSize', 16)
        else
            plot(ax(i), app.(tablename).Data(ts,1), app.(tablename).Data(ts,2), 'o', 'Color', cololist(rem(ts,6)+1,:))
            text(ax(i), app.(tablename).Data(ts,1)+.1, app.(tablename).Data(ts,2), num2str(ts))
        end
    end
    hold(ax(i), 'off');
    xlim(ax(i), [app.T1Edit.Value app.T2Edit.Value])
    ylim(ax(i), [app.Y1Edit.Value-app.Y2Edit.Value app.Y1Edit.Value+app.Y2Edit.Value])
end
app.TSDTable.Data = [(1:12)', app.TS2Table.Data - app.TS1Table.Data];
plot(app.AnalysisAxes, 0, 0);
hold(app.AnalysisAxes, 'on');
plot(app.AnalysisAxes, t, y(:,1), 'b');
plot(app.AnalysisAxes, t, y(:,2), 'r');
for ts = 1:12
    plot(app.AnalysisAxes, app.TS1Table.Data(ts,1), app.TS1Table.Data(ts,2), 'bo')
    text(app.AnalysisAxes, app.TS1Table.Data(ts,1)+.1, app.TS1Table.Data(ts,2), num2str(ts))
    plot(app.AnalysisAxes, app.TS2Table.Data(ts,1), app.TS2Table.Data(ts,2), 'ro')
    text(app.AnalysisAxes, app.TS2Table.Data(ts,1)+.1, app.TS2Table.Data(ts,2), num2str(ts))
%     quiver(app.AnalysisAxes, app.TS1Table.Data(ts,1), app.TS1Table.Data(ts,2), ...
%         app.TS2Table.Data(ts,1)-app.TS1Table.Data(ts,1), app.TS2Table.Data(ts,2)-app.TS1Table.Data(ts,2), 0, ...
%         'k', 'LineWidth', 1, 'MaxHeadSize', 2, 'AutoScale', 'off')
end
hold(app.AnalysisAxes, 'off');
xlim(app.AnalysisAxes, [app.T1Edit.Value app.T2Edit.Value])
ylim(app.AnalysisAxes, [app.Y1Edit.Value-app.Y2Edit.Value app.Y1Edit.Value+app.Y2Edit.Value])
