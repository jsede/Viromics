import React from 'react';
import { AlertTriangle, RefreshCw } from 'lucide-react';

export default class ErrorBoundary extends React.Component {
  constructor(props) {
    super(props);
    this.state = { hasError: false, error: null };
  }

  static getDerivedStateFromError(error) {
    return { hasError: true, error };
  }

  componentDidCatch(error, errorInfo) {
    console.error('Dashboard component crashed:', error, errorInfo);
  }

  render() {
    if (this.state.hasError) {
      return (
        <div className="p-8 glass-panel border-[#EBBDB8] bg-[#FCEEEB] text-center flex flex-col items-center gap-4 my-6">
          <div className="w-10 h-10 rounded-[4px] bg-white border border-[#EBBDB8] flex items-center justify-center">
            <AlertTriangle className="w-5 h-5 text-[#A62C22]" />
          </div>
          <div>
            <h3 className="text-[15px] font-semibold text-[#111827]">Component Display Error</h3>
            <p className="text-[13px] text-[#4B5563] mt-1 max-w-md">
              An error occurred while rendering this tab ({this.props.componentName || 'Active View'}).
            </p>
            {this.state.error && (
              <pre className="mt-3 p-3 bg-white border border-[#E2E6EA] rounded-[4px] text-left text-[11px] font-mono text-[#A62C22] max-w-lg overflow-x-auto select-text">
                {this.state.error.message || String(this.state.error)}
              </pre>
            )}
          </div>
          <button
            onClick={() => this.setState({ hasError: false, error: null })}
            className="flex items-center gap-1.5 px-3.5 py-1.5 bg-white border border-[#CBD2D9] text-xs text-[#111827] font-semibold rounded-[4px] hover:border-[#0F5C6E] hover:text-[#0F5C6E] transition-colors"
          >
            <RefreshCw className="w-3.5 h-3.5 text-[#0F5C6E]" />
            Retry View
          </button>
        </div>
      );
    }
    return this.props.children;
  }
}
