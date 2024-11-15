import ReactMarkdown from "react-markdown";
import remarkGfm from "remark-gfm";
import helpMarkdown from "./Help.md?raw";

const HelpPage = () => {
    return (
        <ReactMarkdown
            className="markdown"
            remarkPlugins={[[remarkGfm]]}
            components={{
                /* 
                react-markdown doesn't support subscript.
                It converts content within ~~ to be within del tags.
                So instead we convert the del tags to sub tags.
                See https://github.com/remarkjs/react-markdown?tab=readme-ov-file#appendix-b-components
                */
                del(props) {
                    const { node, ...rest } = props;
                    return <sub {...rest} />;
                },
            }}
        >
            {helpMarkdown}
        </ReactMarkdown>
    );
};

export default HelpPage;
