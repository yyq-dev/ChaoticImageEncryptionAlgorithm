classdef guiencry < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure       matlab.ui.Figure
        AboutButton    matlab.ui.control.Button
        DecryptButton  matlab.ui.control.Button
        EncryptButton  matlab.ui.control.Button
        SelectButton   matlab.ui.control.Button
        Label          matlab.ui.control.Label
        Label_2        matlab.ui.control.Label
        UIAxes_dec     matlab.ui.control.UIAxes
        UIAxes_enc     matlab.ui.control.UIAxes
        UIAxes_ori     matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        im_gray % Description
        K
        im_ency
        im_decy        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Callback function
        function ImageClicked(app, event)
            
        end

        % Button pushed function: SelectButton
        function SelectButtonPushed(app, event)
            [file,path]=uigetfile({'*.png';'*.jpg';'*.bmp';'*.tiff';'*.*'},'选择图像');
            if isequal(file,0)||isequal(path,0)
                errordlg('未选择图像','错误')
            else
                imfile=strcat(path,file);
            end

            im_orin=imread(imfile);
            app.im_gray=rgb2gray(im_orin);
            imshow(app.im_gray,'Parent',app.UIAxes_ori)
        end

        % Button pushed function: EncryptButton
        function EncryptButtonPushed(app, event)
            options.Resize = 'on';
            options.WindowStyle = 'normal';
            options.FontSize = 12;
            ak=inputdlg('请输入密钥(4行)','请输入密钥',4,{''},options);
            a=cell2mat(ak);
            app.K=str2num(a);
            app.im_ency=tpencrypt(app.im_gray,app.K);
            imshow(app.im_ency,'Parent',app.UIAxes_enc);
        end

        % Button pushed function: DecryptButton
        function DecryptButtonPushed(app, event)
            options.Resize = 'on';
            options.WindowStyle = 'normal';
            options.FontSize = 12;
            ak=inputdlg('请输入密钥(4行)','请输入密钥',4,{''},options);
            a=cell2mat(ak);
            app.K=str2num(a);            
            app.im_decy=tpdecrypt(app.im_ency,app.K);
            imshow(app.im_decy,'Parent',app.UIAxes_dec);
        end

        % Button pushed function: AboutButton
        function AboutButtonPushed(app, event)
            msgbox('版本V1.0.240519','关于','help');
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.9412 0.9412 0.9412];
            app.UIFigure.Position = [100 100 698 445];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes_ori
            app.UIAxes_ori = uiaxes(app.UIFigure);
            title(app.UIAxes_ori, '灰度图')
            app.UIAxes_ori.FontName = '苹方';
            app.UIAxes_ori.XTick = [];
            app.UIAxes_ori.YTick = [];
            app.UIAxes_ori.TitleFontWeight = 'bold';
            app.UIAxes_ori.Position = [1 114 221 218];

            % Create UIAxes_enc
            app.UIAxes_enc = uiaxes(app.UIFigure);
            title(app.UIAxes_enc, '加密后')
            app.UIAxes_enc.FontName = '苹方';
            app.UIAxes_enc.XTick = [];
            app.UIAxes_enc.YTick = [];
            app.UIAxes_enc.TitleFontWeight = 'bold';
            app.UIAxes_enc.Position = [229 114 221 218];

            % Create UIAxes_dec
            app.UIAxes_dec = uiaxes(app.UIFigure);
            title(app.UIAxes_dec, '解密后')
            app.UIAxes_dec.FontName = '苹方';
            app.UIAxes_dec.XTick = [];
            app.UIAxes_dec.YTick = [];
            app.UIAxes_dec.TitleFontWeight = 'bold';
            app.UIAxes_dec.Position = [456 114 221 218];

            % Create Label_2
            app.Label_2 = uilabel(app.UIFigure);
            app.Label_2.BackgroundColor = [0.0745 0.6235 1];
            app.Label_2.HorizontalAlignment = 'center';
            app.Label_2.FontName = '苹方';
            app.Label_2.FontSize = 24;
            app.Label_2.FontWeight = 'bold';
            app.Label_2.FontColor = [1 1 1];
            app.Label_2.Position = [1 1 698 114];
            app.Label_2.Text = '';

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.BackgroundColor = [0.0745 0.6235 1];
            app.Label.HorizontalAlignment = 'center';
            app.Label.FontName = '苹方';
            app.Label.FontSize = 24;
            app.Label.FontWeight = 'bold';
            app.Label.FontColor = [1 1 1];
            app.Label.Position = [1 365 698 81];
            app.Label.Text = '电 力 巡 检 图 像 加 密 系 统';

            % Create SelectButton
            app.SelectButton = uibutton(app.UIFigure, 'push');
            app.SelectButton.ButtonPushedFcn = createCallbackFcn(app, @SelectButtonPushed, true);
            app.SelectButton.FontName = '苹方';
            app.SelectButton.Position = [62 43 100 30];
            app.SelectButton.Text = '打 开';

            % Create EncryptButton
            app.EncryptButton = uibutton(app.UIFigure, 'push');
            app.EncryptButton.ButtonPushedFcn = createCallbackFcn(app, @EncryptButtonPushed, true);
            app.EncryptButton.FontName = '苹方';
            app.EncryptButton.Position = [221 43 100 30];
            app.EncryptButton.Text = '加 密';

            % Create DecryptButton
            app.DecryptButton = uibutton(app.UIFigure, 'push');
            app.DecryptButton.ButtonPushedFcn = createCallbackFcn(app, @DecryptButtonPushed, true);
            app.DecryptButton.FontName = '苹方';
            app.DecryptButton.Position = [384 43 100 30];
            app.DecryptButton.Text = '解 密';

            % Create AboutButton
            app.AboutButton = uibutton(app.UIFigure, 'push');
            app.AboutButton.ButtonPushedFcn = createCallbackFcn(app, @AboutButtonPushed, true);
            app.AboutButton.FontName = '苹方';
            app.AboutButton.Position = [543 43 100 30];
            app.AboutButton.Text = '关 于';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = guiencry

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end