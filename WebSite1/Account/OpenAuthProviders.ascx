<%@ Control Language="C#" AutoEventWireup="true" CodeFile="OpenAuthProviders.ascx.cs" Inherits="OpenAuthProviders" %>

<div id="socialLoginList">
    <h4>Zaloguj się za pomocą innej usługi.</h4>
    <hr />
    <asp:ListView runat="server" ID="providerDetails" ItemType="System.String"
        SelectMethod="GetProviderNames" ViewStateMode="Disabled">
        <ItemTemplate>
            <p>
                <button type="submit" class="btn btn-default" name="provider" value="<%#: Item %>"
                    title="Zaloguj się za pomocą swojego <%#: Item %> konta.">
                    <%#: Item %>
                </button>
            </p>
        </ItemTemplate>
        <EmptyDataTemplate>
            <div>
                <p>Nie skonfigurowano żadnych zewnętrznych usług uwierzytelniania. Zapoznaj się z <a href="https://go.microsoft.com/fwlink/?LinkId=252803">tym artykułem</a>, aby uzyskać szczegółowe informacje dotyczące konfigurowania tej aplikacji platformy ASP.NET w celu obsługi logowania za pośrednictwem usług zewnętrznych.</p>
            </div>
        </EmptyDataTemplate>
    </asp:ListView>
</div>